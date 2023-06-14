import numpy as np
from pyscf import gto, tools, dft
import os, sys
import configparser
from pyscf.geomopt.geometric_solver import optimize

#functions

#saving cube files for homo and lumo orbtials
def save_cubefile(folder_path, homo_index, lumo_index):
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_homo.cube', orbitals[:, homo_index])
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_homo_minus.cube', orbitals[:, homo_index - 1])
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_lumo.cube', orbitals[:, lumo_index])
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_lumo_plus.cube', orbitals[:, lumo_index + 1])

#saving orbitals
def save_orbitals(folder_path, homo_index, lumo_index):
    np.savetxt(folder_path + '/' + 'homo.txt', orbitals[:, homo_index])
    np.savetxt(folder_path + '/' + 'homo_minus.txt', orbitals[:, homo_index - 1])
    np.savetxt(folder_path + '/' + 'lumo.txt', orbitals[:, lumo_index])
    np.savetxt(folder_path + '/' + 'lumo_plus.txt', orbitals[:, homo_index + 1])

#saving optimized geometry as xyz-file
def save_opt_xyz_format(mol, system, folder_path = ''): # added
    coords = mol.atom_coords(unit='ang')
    with open(folder_path + '/' + system + '_optimized.xyz', 'w') as xyz_file:
        xyz_file.write(str(mol.natm) + '\n')
        xyz_file.write(system + '\n')
        j = 0
        for atom in coords:
            list = []
            list.append(mol.atom_symbol(j))
            for i in range(len(atom)):
                list.append(atom[i])
            s = '      '.join(str(x) for x in list)
            j += 1
            xyz_file.write(s+'\n')


if __name__ == "__main__":

    # Read config file

    config = configparser.ConfigParser()
    config.read(sys.argv[2])
    functionals = config['functionals']
    basis_sets = config['basis']
    system = config['system']['molecule'] #added

    # initialize molecule
    mol = gto.Mole()
    mol.atom = sys.argv[1]

    # calculations with all configurations

    for basis in basis_sets:
        for functional in functionals:
            mol.basis = basis_sets[basis]
            mol.build()

            # initialize dft instance
            mf = dft.RKS(mol, xc=functionals[functional])

            # geometry optimized calculation
            mol_opt = optimize(mf)
            mf_opt = dft.RKS(mol_opt, xc=functionals[functional])
            mf_opt.kernel()

            # get orbitals and orbital energies
            energies = mf_opt.mo_energy
            orbitals = mf_opt.mo_coeff

            # get LUMO and HOMO index; ENERGIES ARE ALREADY SORTED!
            lumo_idx = mf_opt.mo_occ.tolist().index(0.0)
            homo_idx = lumo_idx - 1

            # homo-lumo energies and gap
            homo_energy = energies[homo_idx]
            lumo_energy = energies[lumo_idx]
            homo_lumo_gap = lumo_energy - homo_energy

            # SAVING DATA

            # define destination folder
            output_folder = './testsystems/' + system + '/' + functionals[functional] + '_' + \
                            basis_sets[basis] + '/results'

            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

            # save optimized geometry as xyz-file
            save_opt_xyz_format(mol_opt, system, output_folder) #added

            # save energies
            np.savetxt(output_folder + '/' + 'energies.txt', energies)
            homo_lumo_filename = 'homo_lumo_energies.txt'
            with open(output_folder + '/' + homo_lumo_filename, 'w') as file:
                file.write("Total energy: {:.6f} Hartree\n".format(mf_opt.e_tot))
                file.write("HOMO energy: {:.6f} Hartree\n".format(homo_energy)) #check if its really Hartree
                file.write("LUMO energy: {:.6f} Hartree\n".format(lumo_energy))
                file.write("HOMO-LUMO gap: {:.6f} Hartree\n".format(homo_lumo_gap))

            # saving cube files for homo and lumo orbtials
            save_cubefile(output_folder, homo_idx, lumo_idx)
            save_orbitals(output_folder, homo_idx, lumo_idx)


