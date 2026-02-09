#!/usr/bin/env python3

import argparse
import os
import math
import hoomd
import hoomd.md
import numpy as np
import time
def tag(x, n): return f"{x:.{n}f}"
def mic_displacement_ortho(r, r0, L):
        """Minimum-image displacement for orthorhombic boxes.
        r, r0: (N,3) positions; L: (3,) box lengths"""
        dr = r - r0

        dr -= L * np.round(dr / L)
        return dr
class ActiveLJPairs(hoomd.custom.Action):

    def __init__(self, sim, lj):
        self.sim = sim
        self.lj = lj

        self._lj_pairs = 0


        # Build bond exclusion set
        snap = sim.state.get_snapshot()
        self.bond_pairs = set()

        if snap.communicator.rank == 0 and snap.bonds.N > 0:
            for i, j in snap.bonds.group:
                self.bond_pairs.add((min(i, j), max(i, j)))
    
    def act(self, timestep):
        snap = self.sim.state.get_snapshot()
        if snap.communicator.rank != 0:
            return

        pos = snap.particles.position
        types = snap.particles.types
        typeid = snap.particles.typeid
        box = snap.configuration.box  # [Lx, Ly, Lz, xy, xz, yz]
        L = np.array(box[:3], dtype=np.float64)

        N = len(pos)
        count = 0

        for i in range(N):
            ti = types[typeid[i]]
            for j in range(i + 1, N):

                # EXCLUDE bonded pairs
                if (i, j) in self.bond_pairs:
                    continue

                tj = types[typeid[j]]
                rcut = self.lj.r_cut[(ti, tj)]
                if rcut <= 0:
                    continue

                rij = mic_displacement_ortho(np.array(pos[j], dtype=np.float64),np.array(pos[i], dtype=np.float64), L)
                if np.dot(rij, rij) < rcut * rcut:
                    count += 1

        self._lj_pairs = count



    @hoomd.logging.log(category="scalar", default=True)
    def active_lj_pairs(self):
        return self._count

class PressureTrace(hoomd.custom.Action):

    def __init__(self, thermo, volume):
        self.thermo = thermo
        self.volume = volume
        self._trace = 0.0

    def act(self, timestep):
        # pressure_tensor is available after integration starts
        self._trace = (self.thermo.pressure * 3 * self.volume - 2 * self.thermo.kinetic_energy) / self.volume

    @hoomd.logging.log(category="scalar", default=True)
    def pressure_trace(self):
        return self._trace
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
    ap.add_argument("--rho", type=float, default=10.0, help="Number density (segments per volume).")
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
    ap.add_argument("--fluid_dir", type=str, default="/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/snapshotSeed")
    args = ap.parse_args()


    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------ INIT DEVICE ------------------------
    device = hoomd.device.auto_select()
    sim = hoomd.Simulation(device=device, seed=args.seed)
    fluidname = os.path.join(
            args.fluid_dir,
            f"snapshotSeed_rho{tag(args.rho,6)}_NVT.gsd"
        )
    sim.create_state_from_gsd(
        filename=fluidname,
        frame=0   # or -1 for last frame
    )
    box_V = 900 * 25 / args.rho
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

    
    mttk_production = hoomd.md.methods.thermostats.MTTK(
        kT=args.kT,
        tau=args.mttk_tau,
    )

    nvt_production = hoomd.md.methods.ConstantVolume(
        filter=hoomd.filter.All(),
        thermostat=mttk_production
    )

    sim.operations.integrator.methods.append(nvt_production)
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=args.kT)
   
    sim.operations.integrator.forces.append(lj)
    sim.operations.integrator.forces.append(rev_cross)
    sim.run(500000)
    fname = os.path.join(
            args.outdir,
            f"fluidphase_rho{tag(args.rho,6)}_NVT.gsd"
        )
    w = hoomd.write.GSD(
        filename=fname,
        trigger=hoomd.trigger.Periodic(500000),  # write next step only
        mode="ab",
        filter=hoomd.filter.All(),
    )
    sim.operations.writers.append(w)
    sim.run(50000001)
if __name__ == "__main__":
    main()