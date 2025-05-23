# MACE Protein Modelling

## Setup scripts to create solvated systems, systems in vacuum and and assessing NPT equillibration

`aqueous_protein_setup.ipynb`: notebook is where the solvated and vacuum systems were set up for the protein. The initial NPT equillibration was also assessed here.

`solvated_750_water_molecules_2jof_trp_cage.pdb`: solvated system pdb file.

`tip3p.pdb`: TIP3P water system for solvation.

`2jof_trp_cage.pdb`: different conformations from experimental retro Trp-cage structures.

`2luf.pdb`: different conformations from experimental retro Trp-cage structures.

`fully_optimize.py`: python script to optimize initial system setups using LFBGS optimizer.


## Scripts to run Molecular Dynamics simulations

`mace_npt_protein_with_equillibration`: NPT equilibration run.

`continue_run_mace_protein_with_equillibration.py`: regular simulation at constant temperature.

`mace_protein_unfold.py`: simulation at extreme temperatures, increase temperature from starting temperature to ensure system stability.


## Analysis Scripts

`folding_and_unfolding_analysis.ipynb`: main analysis of different trajectories for solvated and vacuum system, calculating radius of gyration and RMSD over the simulation length. Visualizing the trajectory with NGLView to create snapshots for figures.

`energies_and_further_structural_analysis.ipynb`: looking at energies over simulation, producing free energy landscape plot, RMSD vacuum plot, other structural analysis for processing.



