from pyscf import gto
from pyscf import dft
import configparser
import os
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

mf = dft.RKS(methane)
mf.xc = 'pbe'
output_folder = '/Users/severinkeller/Desktop/test.txt'
with open(output_folder, 'w') as f:
    f.write(str(mf.analyze))



