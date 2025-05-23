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

#Include D3 dispersion interactions

macemp = mace_mp(model="small", dispersion=True, device = "cuda", enable_cueq=True)


# Load the PDB structure for the solvated protein
protein = read('/data/fast-pc-02/CRSid/foundation_mace_assignment/vacuum_retro_trp_cage.xyz', index = 0)

print(len(protein))
print(protein)

# Set up ASE calculator
protein.calc = macemp

# Define set up parameters
time_step = 1
T_init = 298
T_target = 498
n_steps_main = 2000000 # See unfolding
friction_equillibration = 0.1
friction_main = 0.01
write_interval = 100

# Initialize velocities according to the Boltzmann distribution

MaxwellBoltzmannDistribution(atoms = protein, temperature_K = T_init)


# Set up the Langevin dynamics engine for NVT ensemble for equillibration steps.

dyn = Langevin(atoms = protein, timestep = (time_step * units.fs), temperature_K = T_init, friction = friction_equillibration)


# Run simulation

print("Starting Simulation")

start = time.time()

# Initial equillibration period with stronger thermostat

for temp in range(T_init, T_target+1, 5):  # Increase by 5K each step
    
    print("Current thermostat temperature in the system: ", temp)
    
    dyn = Langevin(atoms = protein, timestep = (time_step * units.fs), temperature_K = temp, friction = friction_equillibration)
    
    traj = Trajectory('/data/fast-pc-02/CRSid/foundation_mace_assignment/vacuum_unfold/vacuum_nvt_equi_1_simulation_498K_trp_cage_retro.traj', 'a', protein)
    dyn.attach(traj.write, interval = write_interval)
    
    md_logger = MDLogger(dyn, protein, "/data/fast-pc-02/CRSid/foundation_mace_assignment/vacuum_unfold/vacuum_nvt_equi_1_md_log_498K_trp_cage_retro.txt", header=False, mode = "a", stress=True, peratom=True)
    dyn.attach(md_logger, interval = write_interval)

    dyn.run(1000)  

# Main run with weaker thermostat to allow more physical interactions

dyn = Langevin(atoms = protein, timestep = (1 * units.fs), temperature_K = T_target, friction = friction_main)


# Attach a trajectory writer to save trajectory over simulation length
traj_prod = Trajectory('/data/fast-pc-02/CRSid/foundation_mace_assignment/vacuum_unfold/vacuum_nvt_prod_1_simulation_498K_trp_cage_retro.traj', 'w', protein)
dyn.attach(traj_prod.write, interval = 1000)

                                                                                                                                                                                                               # Set up logging for energy and forces
md_logger_prod = MDLogger(dyn, protein, "/data/fast-pc-02/CRSid/foundation_mace_assignment/vacuum_unfold/vacuum_nvt_prod_1_md_log_498K_trp_cage_retro.txt", header=True, stress=True, peratom=True)

dyn.attach(md_logger_prod, interval = 1000)

dyn.run(n_steps_main)

finish = time.time()

print("The time of execution of above program is :",(finish-start), "seconds")
