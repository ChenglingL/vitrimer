import argparse
import os
import math
import hoomd
import hoomd.md
import numpy as np
import time

device = hoomd.device.auto_select()
sim = hoomd.Simulation(device=device, seed=0)

sim.create_state_from_gsd(
        filename=f'/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/GPUfixed/rho1.1/traj_behemoth_T0.00100000_rho1.100000_dt0.001_wait200000.gsd',
        frame=0   # or -1 for last frame
    )

bond = hoomd.md.bond.Harmonic()
bond.params['NN'] = dict(k=1.0, r0=1.0)
sim.operations.integrator = hoomd.md.Integrator(dt=0.001)
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
    "sigma": 0.5, "n": 10, "epsilon": 1, "lambda3": 1}


mttk_production = hoomd.md.methods.thermostats.MTTK(
    kT=1.0,
    tau=1.0,
)

nvt_production = hoomd.md.methods.ConstantVolume(
    filter=hoomd.filter.All(),
    thermostat=mttk_production
)

sim.operations.integrator.methods.append(nvt_production)
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=1.0)

sim.operations.integrator.forces.append(lj)
sim.operations.integrator.forces.append(rev_cross)

sim.run(0)
with sim.operations.integrator.forces[2].cpu_local_force_arrays as arrays:
    F = arrays.force.copy()
    np.savetxt(
        "/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/GPUfixed/rho1.1/test_Rev_forces_GPU.csv",
        F,
        delimiter=",",
        header="fx,fy,fz",
        comments=""
    )
with sim.operations.integrator.forces[1].cpu_local_force_arrays as arrays:
    F = arrays.force.copy()
    np.savetxt(
        "/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/GPUfixed/rho1.1/test_LJ_forces_GPU.csv",
        F,
        delimiter=",",
        header="fx,fy,fz",
        comments=""
    )
with sim.operations.integrator.forces[0].cpu_local_force_arrays as arrays:
    F = arrays.force.copy()
    np.savetxt(
        "/home/cli428/vitrimer/data/test/vitrimerPaper/NVT/V4Test/GPUfixed/rho1.1/test_spring_forces_GPU.csv",
        F,
        delimiter=",",
        header="fx,fy,fz",
        comments=""
    )
