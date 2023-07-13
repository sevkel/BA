import numpy as np
from pyscf import gto, tools, dft
import os
import sys
import configparser
from pyscf.geomopt.geometric_solver import optimize
import density_functional_approximation_dm21 as dm21
from pylibnxc.pyscf import UKS


# FUNCTIONS

def save_cubefile_alpha_beta(mol_obj, mf_obj, folder_path, homo_alpha_index, lumo_alpha_index, homo_beta_index, lumo_beta_index):

    """saving cube files for homo and lumo orbitals"""

    alpha_coeff = mf_obj.mo_coeff[0,:,:]
    beta_coeff = mf_obj.mo_coeff[1,:,:]

    tools.cubegen.orbital(mol_obj, folder_path + '/alpha_orbital_homo.cube', alpha_coeff[:, homo_alpha_index])
    tools.cubegen.orbital(mol_obj, folder_path + '/alpha_orbital_homo_minus.cube', alpha_coeff[:, homo_alpha_index - 1])
    tools.cubegen.orbital(mol_obj, folder_path + '/alpha_orbital_lumo.cube', alpha_coeff[:, lumo_alpha_index])
    tools.cubegen.orbital(mol_obj, folder_path + '/alpha_orbital_lumo_plus.cube', alpha_coeff[:, lumo_alpha_index + 1])

    tools.cubegen.orbital(mol_obj, folder_path + '/beta_orbital_homo.cube', beta_coeff[:, homo_beta_index])
    tools.cubegen.orbital(mol_obj, folder_path + '/beta_orbital_homo_minus.cube', beta_coeff[:, homo_beta_index - 1])
    tools.cubegen.orbital(mol_obj, folder_path + '/beta_orbital_lumo.cube', beta_coeff[:, lumo_beta_index])
    tools.cubegen.orbital(mol_obj, folder_path + '/beta_orbital_lumo_plus.cube', beta_coeff[:, lumo_beta_index + 1])

def save_orbitals(mf_obj, folder_path, homo_alpha_index, lumo_alpha_index, homo_beta_index, lumo_beta_index):

    """saving orbitals"""

    alpha_coeff = mf_obj.mo_coeff[0,:,:]
    beta_coeff = mf_obj.mo_coeff[1,:,:]
    np.savetxt(folder_path + '/' + 'alpha_homo_coeff.txt', alpha_coeff[:, homo_alpha_index])
    np.savetxt(folder_path + '/' + 'alpha_homo_minus_coeff.txt', alpha_coeff[:, homo_alpha_index - 1])
    np.savetxt(folder_path + '/' + 'alpha_lumo_coeff.txt', alpha_coeff[:, lumo_alpha_index])
    np.savetxt(folder_path + '/' + 'alpha_lumo_plus_coeff.txt', alpha_coeff[:, lumo_alpha_index + 1])

    np.savetxt(folder_path + '/' + 'beta_homo_coeff.txt', beta_coeff[:, homo_beta_index])
    np.savetxt(folder_path + '/' + 'beta_homo_minus_coeff.txt', beta_coeff[:, homo_beta_index - 1])
    np.savetxt(folder_path + '/' + 'beta_lumo_coeff.txt', beta_coeff[:, lumo_beta_index])
    np.savetxt(folder_path + '/' + 'beta_lumo_plus_coeff.txt', beta_coeff[:, lumo_beta_index + 1])

def save_opt_xyz_format(mol_obj, system, folder_path=''):

    """saving optimized geometry as xyz-file"""

    coords = mol_obj.atom_coords(unit='ang')
    with open(folder_path + '/' + system + '_optimized.xyz', 'w') as xyz_file:
        xyz_file.write(str(mol_obj.natm) + '\n')
        xyz_file.write(system + '\n')
        j = 0
        for atom in coords:
            list = []
            list.append(mol_obj.atom_symbol(j))
            for i in range(len(atom)):
                list.append(atom[i])
            s = '      '.join(str(x) for x in list)
            j += 1
            xyz_file.write(s + '\n')

def get_atom_energie(atom_symbol, func, basis):  # for now only for W4-11 set

    """get total atom energy"""

    filepath = './results/atoms_W4-11/' + atom_symbol + '/' + func + '_' + basis + '/energies.txt'
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if 'Total energy' in line:
                parts = line.split(':')
                if len(parts) > 1:
                    energy = float(parts[1].split()[0])
                    return energy

def get_react_energy(mf_obj, mol_obj, func, basis):

    """get molecule reaction energy"""

    atom_symbols = [mol_obj.atom_symbol(i).lower() for i in range(mol_obj.natm)]
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
    mf = dft.UKS(mol)

    """CALCULATIONS WITH ALL CONFIGURATIONS"""

    for basis in basis_sets:
        for functional in functionals:

            """define output folder, depending if atom or molecule. We dont need geometry opt. for single atoms
            for now only for W4-11 set"""

            if config.getboolean('system', 'single_atm'):
                output_folder = './results/' + '/atoms_W4-11/' + system + '/' + functionals[functional] + '_' + \
                                basis_sets[basis]
            else:
                output_folder = './results/' + '/molecs_W4-11/' + system + '/' + functionals[functional] + '_' + \
                                basis_sets[basis]
            if os.path.exists(output_folder):
                continue
            else:
                os.makedirs(output_folder)


            if not config.getboolean('system', 'single_atm'):  # geom opt

                if functionals[functional] in ['DM21', 'DM21m', 'DM21mc', 'DM21mu']:  # treat DM21 functional differently
                    mf.xc = 'pbe0'
                    mol_opt = optimize(mf)  # geometry optimization
                    mf_opt = dft.UKS(mol_opt)
                    if functionals[functional] == 'DM21':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
                    elif functionals[functional] == 'DM21m':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21m)
                    elif functionals[functional] == 'DM21mc':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21mc)
                    elif functionals[functional] == 'DM21mu':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21mu)
                elif functionals[functional] in ['GGA_XC_PBE', 'MGGA_XC_SCAN']:
                    mf = UKS(mol, nxc=functionals[functional], nxc_kind='grid')
                    mol_opt = optimize(mf)
                    mf_opt = UKS(mol_opt, nxc=functionals[functional], nxc_kind='grid')
                else:
                    mf.xc = functionals[functional]
                    mol_opt = optimize(mf)
                    mf_opt = dft.UKS(mol_opt)
                    mf_opt.xc = functionals[functional]

                mf_opt.kernel()

                # get orbital energies
                mo_energies = mf_opt.mo_energy

                # get LUMO and HOMO index for alpha and beta; ENERGIES ARE ALREADY SORTED!
                occ_orbs = mf_opt.mo_occ.tolist()

                alpha_lumo_idx = occ_orbs[0].index(0.0)
                beta_lumo_idx = occ_orbs[1].index(0.0)

                alpha_homo_idx = alpha_lumo_idx - 1
                beta_homo_idx = beta_lumo_idx - 1

                # total energy, homo-lumo energies and gap and reaction energy for molecules
                e_tot = mf_opt.e_tot
                energy_alpha_homo = mo_energies[0][alpha_homo_idx]
                energy_alpha_lumo = mo_energies[0][alpha_lumo_idx]
                energy_beta_homo = mo_energies[1][beta_homo_idx]
                energy_beta_lumo = mo_energies[1][beta_lumo_idx]

                alpha_hl_gap = abs(energy_alpha_lumo - energy_alpha_homo)
                beta_hl_gap = abs(energy_beta_lumo - energy_beta_homo)

                e_react = get_react_energy(mf_opt, mol_opt, functionals[functional], basis_sets[basis])

                # save optimized geometry
                save_opt_xyz_format(mol_opt, system, output_folder)

            else:  # no geom opt --> single atoms

                if functionals[functional] == 'DM21':  # treat DM21 functional differently
                    mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
                elif functionals[functional] == 'DM21m':
                    mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21m)
                elif functionals[functional] == 'DM21mc':
                    mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mc)
                elif functionals[functional] == 'DM21mu':
                    mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mu)
                elif functionals[functional] in ['GGA_XC_PBE', 'MGGA_XC_SCAN']:
                    mf = UKS(mol, nxc=functionals[functional], nxc_kind='grid')
                else:
                    mf.xc = functionals[functional]

                mf.kernel()

                # total energy and mo energies
                mo_energies = mf.mo_energy
                e_tot = mf.e_tot


            # SAVING DATA

            # save energies
            np.savetxt(output_folder + '/' + 'mo_energies.txt', mo_energies)
            energy_filename = 'energies.txt'
            with open(output_folder + '/' + energy_filename, 'w') as file:
                file.write("Total energy: {:.15f} Hartree\n".format(e_tot))
                if not config.getboolean('system', 'single_atm'):
                    file.write("Reaction energy: {:.15f} kcal/mol\n".format(e_react * 627.5096080305927))
                    file.write("alpha-HOMO energy: {:.15f} Hartree\n".format(energy_alpha_homo))
                    file.write("alpha-LUMO energy: {:.15f} Hartree\n".format(energy_alpha_lumo))
                    file.write("beta-HOMO energy: {:.15f} Hartree\n".format(energy_beta_homo))
                    file.write("beta-LUMO energy: {:.15f} Hartree\n".format(energy_beta_lumo))
                    file.write("alpha-HOMO-LUMO gap: {:.15f} Hartree\n".format(alpha_hl_gap))
                    file.write("beta-HOMO-LUMO gap: {:.15f} Hartree\n".format(beta_hl_gap))


            # saving cube files for alpha/beta homo and lumo orbitals
            if not config.getboolean('system', 'single_atm'):
                save_cubefile_alpha_beta(mol, mf_opt, output_folder, alpha_homo_idx, alpha_lumo_idx, beta_homo_idx, beta_lumo_idx)
                save_orbitals(mf_opt, output_folder, alpha_homo_idx, alpha_lumo_idx, beta_homo_idx, beta_lumo_idx)


