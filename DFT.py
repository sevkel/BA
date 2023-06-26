import numpy as np
from pyscf import gto, tools, dft
import os
import sys
import configparser
from pyscf.geomopt.geometric_solver import optimize
import density_functional_approximation_dm21 as dm21


# Functions

# saving cube files for homo and lumo orbitals
def save_cubefile(folder_path, homo_index, lumo_index):
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_homo.cube', orbitals[:, homo_index])
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_homo_minus.cube', orbitals[:, homo_index - 1])
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_lumo.cube', orbitals[:, lumo_index])
    tools.cubegen.orbital(mol, folder_path + '/' + 'orbital_lumo_plus.cube', orbitals[:, lumo_index + 1])


# saving orbitals
def save_orbitals(folder_path, homo_index, lumo_index):
    np.savetxt(folder_path + '/' + 'homo.txt', orbitals[:, homo_index])
    np.savetxt(folder_path + '/' + 'homo_minus.txt', orbitals[:, homo_index - 1])
    np.savetxt(folder_path + '/' + 'lumo.txt', orbitals[:, lumo_index])
    np.savetxt(folder_path + '/' + 'lumo_plus.txt', orbitals[:, homo_index + 1])


# saving optimized geometry as xyz-file
def save_opt_xyz_format(mol, system, folder_path=''):
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
            xyz_file.write(s + '\n')


# get total atom energy

def get_atom_energie(atom_symbol, func, basis):  # for now only for W4-11 set
    filepath = './results/atoms_W4-11/' + atom_symbol + '/' + func + '_' + basis + '/energies.txt'
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if 'Total energy' in line:
                parts = line.split(':')
                if len(parts) > 1:
                    energy = float(parts[1].split()[0])
                    return energy


# get molecule reaction energy

def get_react_energy(mf_obj, mol, func, basis):
    atom_symbols = [mol.atom_symbol(i).lower() for i in range(mol.natm)]
    e_react = -mf_obj.e_tot
    counted_atoms = set()
    for atom in atom_symbols:
        if atom not in counted_atoms:
            occurrence = atom_symbols.count(atom)
            atom_etot = get_atom_energie(atom, func, basis)
            e_react += occurrence * atom_etot
            counted_atoms.add(atom)
    return e_react


if __name__ == "__main__":

    # Read config file

    config = configparser.ConfigParser()
    config.read(sys.argv[2])
    functionals = config['functionals']
    basis_sets = config['basis']
    system = config['system']['molecule']

    # initialize molecule

    mol = gto.Mole()
    mol.atom = sys.argv[1]
    # mol.verbose = 4 # shows details in building the molecule. un-commend if wanted
    mol.basis = basis_sets['basis1']  # define the basis before the loop as we only use def2-qzvp
    mol.spin = int(config['system']['unpaired_elecs'])
    mol.build()
    if mol.spin == 0:  # restricted if no unpaired elecs
        mf = dft.RKS(mol)
    else:
        mf = dft.UKS(mol)  # unrestricted if there are unpaired elecs

    if not config.getboolean('system', 'single_atm'):
        mf.xc = 'pbe0'  # optimize with the more expensive pbe0 functional (not possible to optimize with the DM21)
        mol_opt = optimize(mf)

    # CALCULATIONS WITH ALL CONFIGURATIONS

    for basis in basis_sets:
        for functional in functionals:

            # define output folder, depending if atom or molecule. We dont need geometry opt. for single atoms
            # for now only for W4-11 set
            # maybe get output folder as argv argument??
            if config.getboolean('system', 'single_atm'):
                output_folder = './results/' + '/atoms_W4-11/' + system + '/' + functionals[functional] + '_' + \
                                basis_sets[basis]
            else:
                output_folder = './results/' + '/molecs_W4-11/' + system + '/' + functionals[functional] + '_' + \
                                basis_sets[basis]
            if os.path.exists(output_folder):
                break
            else:
                os.makedirs(output_folder)

            # geometry optimized calculation, if wanted
            if not config.getboolean('system', 'single_atm'):  # geom opt
                #mol_opt = optimize(mf)
                if mol.spin == 0:
                    mf_opt = dft.RKS(mol_opt)
                else:
                    mf_opt = dft.UKS(mol_opt)

                if functionals[functional] == 'DM21':  # added: treat DM21 functional differently
                    mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
                else:
                    mf_opt.xc = functionals[functional]

                mf_opt.kernel()

                # get orbitals and orbital energies
                mo_energies = mf_opt.mo_energy
                orbitals = mf_opt.mo_coeff

                # get LUMO and HOMO index; ENERGIES ARE ALREADY SORTED!
                lumo_idx = mf_opt.mo_occ.tolist().index(0.0)
                homo_idx = lumo_idx - 1

                # total energy, homo-lumo energies and gap and reaction energy for molecules
                e_tot = mf_opt.e_tot
                homo_energy = mo_energies[homo_idx]
                lumo_energy = mo_energies[lumo_idx]
                homo_lumo_gap = lumo_energy - homo_energy
                e_react = get_react_energy(mf_opt, mol_opt, functionals[functional], basis_sets[basis])

                # save optimized geometry
                save_opt_xyz_format(mol_opt, system, output_folder)

            else:  # (no geom opt)
                mf.xc = functionals[functional]
                mf.kernel()

                # get orbitals and orbital energies
                mo_energies = mf.mo_energy
                orbitals = mf.mo_coeff

                # total energy
                e_tot = mf.e_tot
                '''# get LUMO and HOMO index; ENERGIES ARE ALREADY SORTED!
                lumo_idx = mf_opt.mo_occ.tolist().index(0.0) #Gives a ALPHA/BETA list?
                homo_idx = lumo_idx - 1 

                # total energy, homo-lumo energies and gap and reaction energy for molecules
                e_tot = mf_opt.e_tot
                homo_energy = mo_energies[homo_idx]
                lumo_energy = mo_energies[lumo_idx]
                homo_lumo_gap = lumo_energy - homo_energy
                '''



            # SAVING DATA

            # save energies
            np.savetxt(output_folder + '/' + 'mo_energies.txt', mo_energies)
            energy_filename = 'energies.txt'
            with open(output_folder + '/' + energy_filename, 'w') as file:
                file.write("Total energy: {:.15f} Hartree\n".format(e_tot))
                if not config.getboolean('system', 'single_atm'):
                    file.write("Reaction energy: {:.15f} Hartree\n".format(e_react))
                    file.write("HOMO energy: {:.15f} Hartree\n".format(homo_energy))
                    file.write("LUMO energy: {:.15f} Hartree\n".format(lumo_energy))
                    file.write("HOMO-LUMO gap: {:.15f} Hartree\n".format(homo_lumo_gap))


            # saving cube files for homo and lumo orbitals
            if not config.getboolean('system', 'single_atm'):
                save_cubefile(output_folder, homo_idx, lumo_idx)
                save_orbitals(output_folder, homo_idx, lumo_idx)

