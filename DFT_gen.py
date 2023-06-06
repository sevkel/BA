import numpy as np
from pyscf import gto, tools, dft
import sys
import configparser
import os

#try input folder
#input_cfg = './testsystems/' + sys.argv[2]

#read config file
config = configparser.ConfigParser()
config.read(sys.argv[2])
functional = config['DFT']['functional']
basis = config['DFT']['basis']

#initialize molecule
mol = gto.Mole()
mol.atom = sys.argv[1]
mol.basis = basis
mol.build()

#DFT calculation with XC-functional
mf = dft.RKS(mol)
mf.xc = functional
mf.kernel()

# get orbitals and orbital energies
energies = mf.mo_energy
orbitals = mf.mo_coeff

# get LUMO and HOMO index; ENERGIES ARE ALREADY SORTED!
lumo_idx = mf.mo_occ.tolist().index(0.0)
homo_idx = lumo_idx - 1


#SAVING DATA

output_folder = './testsystems/' + config['SYSTEM']['molecule'] + '/' + functional + '_' + basis + '/results'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#save energies
np.savetxt(output_folder + '/' + 'energies.txt', energies)

#cube files for homo and lumo orbtials
tools.cubegen.orbital(mol, output_folder + '/' + 'orbital_homo.cube', orbitals[:, homo_idx])
tools.cubegen.orbital(mol, output_folder + '/' + 'orbital_homo_minus.cube', orbitals[:, homo_idx-1])
tools.cubegen.orbital(mol, output_folder + '/' + 'orbital_lumo.cube', orbitals[:, lumo_idx])
tools.cubegen.orbital(mol, output_folder + '/' + 'orbital_lumo_plus.cube', orbitals[:, lumo_idx+1])

#save orbitals
np.savetxt(output_folder + '/' + 'homo.txt', orbitals[:, homo_idx])
np.savetxt(output_folder + '/' + 'homo_minus.txt', orbitals[:, homo_idx-1])
np.savetxt(output_folder + '/' + 'lumo.txt', orbitals[:, lumo_idx])
np.savetxt(output_folder + '/' + 'lumo_plus.txt', orbitals[:, homo_idx+1])


