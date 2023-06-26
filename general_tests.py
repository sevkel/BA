import density_functional_approximation_dm21 as dm21
from pyscf import gto
from pyscf import dft
from pyscf.geomopt.geometric_solver import optimize

# Create the molecule of interest and select the basis set.
methane = gto.Mole()
methane.atom = """H 0.000000000000 0.000000000000 0.000000000000
                  C 0.000000000000 0.000000000000 1.087900000000
                  H 1.025681956337 0.000000000000 1.450533333333
                  H -0.512840978169 0.888266630391 1.450533333333
                  H -0.512840978169 -0.888266630391 1.450533333333"""
methane.basis = 'def2-qzvp'
#methane.verbose = 4
methane.build()
'''atom_symbols = [methane.atom_symbol(i).lower() for i in range(methane.natm)]
print(atom_symbols)
counted_atoms = set()
for atom in atom_symbols:
  if atom not in counted_atoms:
    occurrence = atom_symbols.count(atom)
    print(str(occurrence) + ' ' + atom)
    counted_atoms.add(atom)'''

carbon = gto.Mole()
carbon.atom = 'C 0.0 0.0 0.0'
carbon.basis = 'def2-qzvp'
carbon.spin = 2
#carbon.verbose = 4
carbon.build()

hydrogen = gto.Mole()
hydrogen.atom = 'H 0.0 0.0 0.0'
hydrogen.basis = 'def2-qzvp'
hydrogen.spin = 1
#hydrogen.verbose = 4
hydrogen.build()

energies = []

mf_meth = dft.RKS(methane)
meth_opt = optimize(mf_meth)
mf_meth_opt = dft.RKS(meth_opt)
mf_meth_opt.xc = 'B3LYP'
mf_meth_opt.run()
dm0 = mf_meth_opt.make_rdm1()

mf_meth_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
# It's wise to relax convergence tolerances.
mf_meth_opt.conv_tol = 1E-6
mf_meth_opt.conv_tol_grad = 1E-3
# Run the DFT calculation.
energy = mf_meth_opt.kernel(dm0=dm0)
energies.append(energy)

for mol in [carbon, hydrogen]:
  # Create a DFT solver and insert the DM21 functional into the solver.
  if mol.spin == 0:
    mf = dft.RKS(mol)
  else:
    mf = dft.UKS(mol)
  # It will make SCF faster to start close to the solution with a cheaper
  # functional.
  mf.xc = 'B3LYP'
  mf.run()
  dm0 = mf.make_rdm1()

  mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
  # It's wise to relax convergence tolerances.
  mf.conv_tol = 1E-6
  mf.conv_tol_grad = 1E-3
  # Run the DFT calculation.
  energy = mf.kernel(dm0=dm0)
  energies.append(energy)

print({'CH4': energies[0], 'C': energies[1], 'H': energies[2]})