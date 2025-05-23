## Imports

import time

# Mace Libraries

from mace.calculators import mace_mp

# Molecular Dynamics Libraries

from ase import build
from ase.md import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.io import Trajectory
from ase.io import read, write
from ase.visualize import view
from ase.md import MDLogger

# Main Code

#Include D3 dispersion interactions

macemp = mace_mp(model="small", dispersion=True, device = "cuda", enable_cueq=True)


# Load the Protein structure from end point of previous run
previous_run_traj = read("/home/CRSid/rds/hpc-work/written_assignment_1/results/solvent_npt_trajectory.traj", index = ":")

protein = previous_run_traj[-1]

# Set up ASE calculator
protein.calc = macemp


# Define set up parameters
time_step = 1
T_init = 300
n_steps_main = 435000
friction_main = 0.01
write_interval = 1000

# Set up the Langevin dynamics engine for NVT ensemble.
dyn = Langevin(atoms = protein, timestep = (time_step * units.fs), temperature_K = T_init, friction = friction_main)


# Attach a trajectory writer to save trajectory over simulation length
traj = Trajectory(filename = '/home/CRSid/rds/hpc-work/written_assignment_1/results/solvent_npt_trajectory.traj', mode = 'a', atoms = protein)
dyn.attach(traj.write, interval = write_interval)


# Set up logging for energy and forces
md_logger = MDLogger(dyn, protein, "/home/CRSid/rds/hpc-work/written_assignment_1/results/nvt_production_1_md_log_300K_trp_cage_standard.txt", header=False, mode = "a", stress=True, peratom=True)
dyn.attach(md_logger, interval = write_interval)

# Run simulation

print("Starting Simulation")

start = time.time()

dyn.run(n_steps_main)

finish = time.time()

print("The time of execution of above program is :",(finish-start), "seconds")

