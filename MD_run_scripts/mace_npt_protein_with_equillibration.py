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

from ase.md.nptberendsen import NPTBerendsen

# Main Code

#Include D3 dispersion interactions

macemp = mace_mp(model="small", dispersion=True, device = "cuda", enable_cueq=True)


# Load the PDB structure for the solvated protein
protein = read('/home/CRSid/rds/hpc-work/written_assignment_1/npt_start_optimized_centred_0_1_solvated_750_water_molecules_2jof_trp_cage.pdb', index = 0)

print(len(protein))
print(protein)

# Set up ASE calculator
protein.calc = macemp

# Define set up parameters
time_step = 1
T_init = 300
n_steps_equillibration = 5000
friction_equillibration = 0.1
write_interval = 10

# Initialize velocities according to the Boltzmann Distribution.

MaxwellBoltzmannDistribution(atoms = protein, temperature_K = T_init)

# Set up the Berendsen dynamics engine for NPT ensemble.

dyn = NPTBerendsen(protein, timestep = time_step * units.fs, temperature_K = T_init,
                   pressure_au = 1.01325 * units.bar, compressibility_au=4.57e-5 / units.bar)

# Attach a trajectory writer to save trajectory over simulation length
traj = Trajectory('/home/CRSid/rds/hpc-work/written_assignment_1/results/npt_equillibration_atm_300K_trp_cage_standard.traj', 'w', protein)
dyn.attach(traj.write, interval = write_interval)


# Set up logging for energy and forces
md_logger = MDLogger(dyn, protein, "/home/CRSid/rds/hpc-work/written_assignment_1/results/npt_equillibrationmd_log_atm_300K_trp_cage_standard.txt", header=True, stress=True, peratom=True)
dyn.attach(md_logger, interval = write_interval)

# Run simulation

print("Starting Simulation")

start = time.time()

# Initial equillibration period with stronger thermostat

dyn.run(n_steps_equillibration)

finish = time.time()

print("The time of execution of above program is :",(finish-start), "seconds")

