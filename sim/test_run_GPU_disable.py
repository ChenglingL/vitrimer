#!/usr/bin/env python3
## Here we diable different potentials to test different version of v5 and v4
import argparse
import os
import math
import hoomd
import hoomd.md
import numpy as np
import time

# ------------------------ GEOMETRY ------------------------
def random_pos(box_len):
    return np.random.uniform(-box_len / 2, box_len / 2, size=3)

def wrap_into_box(pos, box_lengths):
    """
    Wrap a particle position into the periodic box.

    Parameters
    ----------
    pos : array-like, shape (3,)
        The [x, y, z] position in Cartesian coordinates.
    box_lengths : array-like, shape (3,)
        The box dimensions [Lx, Ly, Lz].

    Returns
    -------
    wrapped_pos : ndarray, shape (3,)
        Position wrapped into the box, in range [-L/2, L/2).
    """
    pos = np.array(pos, dtype=float)
    L = np.array(box_lengths, dtype=float)
    wrapped_pos = ((pos + 0.5 * L) % L) - 0.5 * L
    return wrapped_pos

def create_star(center, arm_types, spacing=1.0, start_index=0):
    positions = [center]
    types = ['C']
    bonds = []
    
    central_idx = start_index  # This will be index 0 in this star
    current_idx = start_index + 1  # Next available index after center

    # Generate 8 directions roughly evenly spaced on a sphere
    golden_angle = np.pi * (3 - np.sqrt(5))
    directions = []
    for i in range(8):
        y = 1 - (i / 7.0) * 2  # y from 1 to -1
        radius = np.sqrt(1 - y * y)
        theta = golden_angle * i
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        directions.append(np.array([x, y, z]))

    for i, vec in enumerate(directions):
        prev = central_idx  # start from the center

        for j in range(3):  # 3 beads per arm
            pos = center + (j + 1) * spacing * vec
            positions.append(pos)

            bond = [prev, current_idx]
            bonds.append(bond)

            if j == 2:
                types.append(arm_types[i])  # Final bead is A or B
            else:
                types.append('M')

            prev = current_idx
            current_idx += 1

    return positions, types, bonds

def base_log_steps(duration, per_decade=20, include_end=True):
    """
    Return integer steps in [1, duration] spaced by Δlog10 = 1/per_decade.
    Example: per_decade=20 -> 20 intervals per decade (0.05 in log10).
    """
    if duration < 1:
        return np.array([], dtype=int)

    # quantize log10 to exact multiples of 1/per_decade
    qmin = 0                       # log10(1) = 0
    qmax = int(np.floor(np.log10(duration) * per_decade))
    if include_end:
        qmax += 1                  # include the upper boundary gridpoint

    grid = (np.arange(qmin, qmax) / per_decade)
    steps = np.rint(10.0**grid).astype(int)
    steps = steps[(steps >= 100) & (steps <= duration)]
    return np.unique(steps)

def steps_for_file(wait, base_steps, end=None):
    """
    Shift base_steps by 'wait'. If end is given, keep steps <= end.
    """
    s = wait + base_steps
    if end is not None:
        s = s[s <= end]
    return s


def write_logspaced_multiple(sim, files, base_steps):
    """
    files: list of dicts with {"filename":..., "wait":..., "end": optional}
    base_steps: 1D array of positive ints (the same for all files)
    """
    # Build step -> list of filenames map
    step_to_files = {}
    for f in files:
        fname = f["filename"]
        wait  = int(f["wait"])
        end   = f.get("end", None)
        file_steps = steps_for_file(wait+sim.timestep, base_steps, end=end)
        for s in file_steps:
            step_to_files.setdefault(int(s), []).append(fname)

    # Run the sim and write exactly at those steps
    current = int(sim.timestep)
    for target in sorted(step_to_files.keys()):
        if target > current:
            sim.run(target - current)
            current = target
        else: continue

        writers = []
        for fname in step_to_files[target]:
            w = hoomd.write.GSD(
                filename=fname,
                trigger=hoomd.trigger.Periodic(1),  # write next step only
                mode="ab",
                filter=hoomd.filter.All(),
            )
            sim.operations.writers.append(w)
            writers.append(w)

        sim.run(1)  # write
        current += 1

        for w in writers:
            sim.operations.writers.remove(w)


# ------------------------ SYSTEM BUILD ------------------------

def main():
    
    print("GPU: \n")
    print(hoomd.device.GPU.get_available_devices())
    ap = argparse.ArgumentParser(description="Star vitrimer (HOOMD v5) runner.")
    ap.add_argument("--kT", type=float, required=True, help="Temperature (kT).")
    ap.add_argument("--rho", type=float, required=True, help="Number density (segments per volume).")
    ap.add_argument("--dt", type=float, default=1e-3)
    ap.add_argument("--mttk_tau", type=float, default=1.0)
    ap.add_argument("--end", type=float, default=None)
    # ap.add_argument("--ep", type=float, default=100.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--duration_after_wait", type=int, default=1000000)
    ap.add_argument("--prod_steps", type=int, default=1_000_000)
    ap.add_argument("--warmup_steps", type=int, default=50_000, help="Soft WCA warm-up steps.")
    ap.add_argument("--per_decade", type=int, default=20)
    ap.add_argument("--waits", type=int, nargs="*", default=[0, 100_000, 300_000],
                    help="List of wait steps for separate files.")
    ap.add_argument("--outdir", type=str, default="/home/cli428/vitrimer/data/test/vitrimerPaper/NVT")
    args = ap.parse_args()

    n_stars_N1 = 600  # 7 A ends, 1 B end
    n_stars_N2 = 300  # 1 A end, 7 B ends
    segments_per_star = 25
    sigma = 0.9
    total_stars = n_stars_N1 + n_stars_N2
    total_segments = total_stars * segments_per_star

    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------ INIT DEVICE ------------------------
    device = hoomd.device.auto_select()
    sim = hoomd.Simulation(device=device, seed=args.seed)
    ini_box_length = 300
    box_vol = total_segments / args.rho
    box_length = (total_segments / args.rho) ** (1 / 3)
    all_pos, all_types, all_bonds = [], [], []
    index_offset = 0

    for _ in range(n_stars_N1):
        arms = ['A'] * 7 + ['B']
        np.random.shuffle(arms)
        pos, types, bonds = create_star(random_pos(ini_box_length), arms,spacing=0.5, start_index=index_offset)
        all_pos.extend(pos)
        all_types.extend(types)
        all_bonds.extend(bonds)
        index_offset += len(pos)

    for _ in range(n_stars_N2):
        arms = ['B'] * 7 + ['A']
        np.random.shuffle(arms)
        pos, types, bonds = create_star(random_pos(ini_box_length), arms,spacing=0.5, start_index=index_offset)
        all_pos.extend(pos)
        all_types.extend(types)
        all_bonds.extend(bonds)
        index_offset += len(pos)

    for i in range(len(all_pos)):
        all_pos[i] = wrap_into_box(all_pos[i],[ini_box_length, ini_box_length, ini_box_length])

    #print("position of [0,0]:", [i for i, pair in enumerate(all_bonds) if pair == [0, 0]])


    # ------------------------ SNAPSHOT ------------------------
    snapshot = hoomd.Snapshot()
    snapshot.particles.N = len(all_pos)
    snapshot.particles.position[:] = all_pos
    snapshot.particles.types = ['C', 'M', 'A', 'B']
    snapshot.particles.typeid[:] = [snapshot.particles.types.index(t) for t in all_types]
    
    snapshot.bonds.N = len(all_bonds)
    snapshot.bonds.types = ['NN']
    snapshot.bonds.group[:] = all_bonds
    snapshot.bonds.typeid[:] = [0] * len(all_bonds)

    volume_ramp = hoomd.variant.box.InverseVolumeRamp([ini_box_length, ini_box_length, ini_box_length, 0, 0, 0], box_vol, 0, 10000)
    snapshot.configuration.box = [ini_box_length, ini_box_length, ini_box_length, 0, 0, 0]
    sim.create_state_from_snapshot(snapshot)

    # ------------------------ BONDED FORCE ------------------------
    bond = hoomd.md.bond.Harmonic()
    bond.params['NN'] = dict(k=1000.0, r0=1.0)
    sim.operations.integrator = hoomd.md.Integrator(dt=args.dt)
    sim.operations.integrator.forces.append(bond)

    # ------------------------ WCA PAIR FORCE ------------------------
    nl = hoomd.md.nlist.Cell(buffer=0.4)
    nl.exclusions = ["bond"]
    lj = hoomd.md.pair.LJ(nlist=nl)
    lj.mode = 'shift'

    wca_cut = 2 ** (1 / 6) * 0.9
    type_names = ['C', 'M', 'A', 'B']
    for t1 in type_names:
        for t2 in type_names:
            if set((t1, t2)) != set(('A', 'B')):
                lj.params[(t1, t2)] = dict(epsilon=1.0, sigma=0.9)
                lj.r_cut[(t1, t2)] = wca_cut
            else:
                lj.params[(t1, t2)] = dict(epsilon=0.0, sigma=0.01)
                lj.r_cut[(t1, t2)] = 0.0
                # lj.r_on[(t1, t2)] = wca_cut + 0.5
                
    #sim.operations.integrator.forces.append(lj)

    # ------------------------ REVERSIBLE CROSSLINKING ------------------------
    rev_cross = hoomd.md.many_body.RevCross(default_r_cut=1.8,nlist=nl)
    for t1 in type_names:
        for t2 in type_names:
            if set((t1, t2)) != set(('A', 'B')):
                rev_cross.params[(t1, t2)] = {"sigma":0,"n": 0, "epsilon": 0, "lambda3": 0}
    rev_cross.params[('A','B')] = {
        "sigma": 0.5, "n": 10, "epsilon": 100, "lambda3": 1}
    # rev_cross.lambda_ = 1.0  # swap barrier control
    # rev_cross.epsilon = 100.0
    #sim.operations.integrator.forces.append(rev_cross)
    #lj.r_cut[('A', 'B')] = 0  # disable A-B in LJ
    #sim.operations.integrator.forces.append(lj)

    # # ------------------------ Soft repulsive force to initialize ------------------------

    gauss = hoomd.md.pair.Gaussian(nlist=nl)

    # Example: all particles interact softly at the start
    for a in ['C', 'M', 'A', 'B']:
        for b in ['C', 'M', 'A', 'B']:
            gauss.params[(a,b)] = dict(epsilon=10.0, sigma=0.4)
            gauss.r_cut[(a,b)] = 3*0.8

    sim.operations.integrator.forces.append(gauss)


        # ------------------------ NVT THERMOSTAT ------------------------
    print("box resize with soft start \n")
    soft_T = 0.2
    langevin = hoomd.md.methods.Langevin(
        kT=soft_T,
        filter=hoomd.filter.All()
    )

    
    
    
    sim.operations.integrator.methods.append(langevin)
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=soft_T)
    
    box_resize = hoomd.update.BoxResize(
        trigger=hoomd.trigger.Periodic(10),  # Example: triggers every 10 timesteps
        box=volume_ramp,  # Pass the variant to define the box change over time
        filter=hoomd.filter.All()
    )

    # 4. Add to simulation
    sim.operations.updaters.append(box_resize)

    sim.run(100000) #soft start
    #sim.run(100)
    print(f"Now box: {sim.state.box} \n")
    sim.operations.updaters.remove(box_resize)
    
    mttk_production = hoomd.md.methods.thermostats.MTTK(
        kT=args.kT,
        tau=args.mttk_tau,
    )

    nvt_production = hoomd.md.methods.ConstantVolume(
        filter=hoomd.filter.All(),
        thermostat=mttk_production
    )
    sim.operations.integrator.methods.remove(langevin)
    sim.operations.integrator.methods.append(nvt_production)
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=args.kT)
   
    sim.operations.integrator.forces.remove(gauss)
    sim.operations.integrator.forces.append(lj)
    #sim.operations.integrator.forces.append(rev_cross)

    def tag(x, n): return f"{x:.{n}f}"
    def mkname(wait):
        return os.path.join(
                args.outdir,
                f"traj_T{tag(args.kT,8)}_rho{tag(args.rho,6)}_dt{tag(args.dt,3)}_wait{wait}.gsd"
        )
    logger = hoomd.logging.Logger(categories=['scalar'])
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)
    logger.add(thermo, quantities=[
        "potential_energy",
        "kinetic_energy",
        "kinetic_temperature"
    ])
    
    f_table = open(os.path.join(
        args.outdir,
        f"energy_T{tag(args.kT,8)}_rho{tag(args.rho,6)}_NVT.log"
    ), "w")
    
    table = hoomd.write.Table(
        logger=logger,
        trigger=hoomd.trigger.Periodic(10000),
        output=f_table
    )
    sim.operations.writers.append(table)

    files = [{"filename": mkname(w), "wait": int(w), "end": args.end} for w in args.waits]
    base = base_log_steps(args.duration_after_wait, per_decade=args.per_decade, include_end=True)
    start_step = sim.timestep
    t0 = time.perf_counter()
    write_logspaced_multiple(sim, files, base)
    t1 = time.perf_counter()
    steps = sim.timestep - start_step
    logfile = open(os.path.join(
        args.outdir,f"timeInfo_T{tag(args.kT,8)}_rho{tag(args.rho,6)}_NVT.log"), "w")
    logfile.write("GPU: ")
    logfile.write(f" {hoomd.device.GPU.get_available_devices()}/n")
    logfile.write(f"steps {steps}\n")
    logfile.write(f"time_sec {t1 - t0:.6f}\n")
    logfile.write(f"steps_per_sec {steps/(t1-t0):.2f}\n")
    logfile.close()
if __name__ == "__main__":
    main()