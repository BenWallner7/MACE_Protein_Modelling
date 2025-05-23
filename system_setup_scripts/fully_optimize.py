from ase.optimize import LBFGS

from ase.io import read, write

from mace.calculators import mace_mp


macemp = mace_mp()

# Load im

atoms = read("protein.pdb")

atoms.calc = macemp

# Perform L-BFGS optimization
optimizer = LBFGS(atoms)

optimizer.run(fmax=0.01)
