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
    mol.basis = basis_sets['basis1']
    mol.spin = int(config['system']['unpaired_elecs'])
    mol.charge = int(config['system']['charge'])
    mol.build()
    mf = dft.UKS(mol).density_fit()



    params_pbe = {
        'convergence_energy': 1e-6,  # Eh
        'convergence_grms': 5e-3,  # Eh/Bohr
        'convergence_gmax': 5e-3,  # Eh/Bohr
        'convergence_drms': 1.2e-2,  # Angstrom
        'convergence_dmax': 1.8e-2,  # Angstrom
    }

    conv_params = {
        'convergence_energy': 1e-6,  # Eh
        'convergence_grms': 3e-4,  # Eh/Bohr
        'convergence_gmax': 4.5e-4,  # Eh/Bohr
        'convergence_drms': 1.2e-3,  # Angstrom
        'convergence_dmax': 1.8e-3,  # Angstrom
    }


    conv_params_rel = {  # relaxed settings
        'convergence_energy': 1e-6,  # Eh
        'convergence_grms': 3e-3,    # Eh/Bohr
        'convergence_gmax': 4.5e-3,  # Eh/Bohr
        'convergence_drms': 2.5e-3,  # Angstrom
        'convergence_dmax': 2.5e-3,  # Angstrom
    }

    """CALCULATIONS WITH ALL CONFIGURATIONS"""

    optimization_flag = True

    for basis in basis_sets:
        for functional in functionals:

            """
            define output folder, depending if atom or molecule. We dont need geometry opt. for single atoms
            """

            if config.getboolean('system', 'single_atm'):
                output_folder = './resultsDFT/' + '/atoms_YOURDATASET/' + system + '/' + functionals[functional] + '_' + \
                                basis_sets[basis]
            else:
                output_folder = './resultsDFT/' + 'molecs_YOURDATASET/' + system + '/' + functionals[functional] + '_' + \
                                basis_sets[basis]
            if os.path.exists(output_folder) and len(os.listdir(output_folder)) > 0:
                continue
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

            if ((config.getboolean('system', 'single_atm') == False) and (optimization_flag == True)):
                mf.xc = 'pbe0'
                mol_pbe0_opt = optimize(mf, **conv_params)
                mf_DM21_pre = dft.UKS(mol_pbe0_opt, xc='B3LYP') # to start closer to the solution with a 'cheaper' functional
                mf_DM21_pre.run()
                dm0 = mf_DM21_pre.make_rdm1()
                optimization_flag = False

            if not config.getboolean('system', 'single_atm'):  # geometry optimization

                if functionals[functional] in ['DM21', 'DM21m', 'DM21mc', 'DM21mu']:
                    mol_opt = mol_pbe0_opt
                    mf_opt = dft.UKS(mol_opt)
                    mf_opt.verbose=4
                    if functionals[functional] == 'DM21':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
                        mf_opt.conv_tol = 1E-6
                        mf_opt.conv_tol_grad = 1E-3
                    elif functionals[functional] == 'DM21m':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21m)
                        mf_opt.conv_tol = 1E-6
                        mf_opt.conv_tol_grad = 1E-3
                    elif functionals[functional] == 'DM21mc':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21mc)
                        mf_opt.conv_tol = 1E-6
                        mf_opt.conv_tol_grad = 1E-3
                    elif functionals[functional] == 'DM21mu':
                        mf_opt._numint = dm21.NeuralNumInt(dm21.Functional.DM21mu)
                        mf_opt.conv_tol = 1E-6
                        mf_opt.conv_tol_grad = 1E-3
                    # mf_opt.damp = 0.5 # uncomment if you reach convergence problems
                    mf_opt.max_cycle = 250
                    mf_opt.kernel(dm0=dm0)
                if functionals[functional] in ['GGA_XC_PBE', 'MGGA_XC_SCAN', 'MGGA_HM', 'GGA_HM']:
                    """
                    use one of the codes below and comment the other, wheter you run into convergence problems with 
                    the functionals or not. 
                    
                    !!! 
                    Also, run this code with only one of these four functionals together 
                    with the other functionals.
                    !!!


                    ### PBE0 relaxed startgeometry for SCF calculation
                    -------------------------
                    mol_opt = mol_pbe0_opt
                    mf_opt = UKS(mol_opt, nxc=functionals[functional], nxc_kind='grid').density_fit()
                    # mf_opt.damp = 0.5 # uncomment if you reach convergence problems
                    mf_opt.max_cycle = 250
                    mf_opt.kernel()
                    -------------------------


                    ### Try relaxation with the corresponding ML-functional
                    -------------------------
                    mf.xc = 'pbe'
                    mf.conv_tol = 1E-6
                    mf.conv_tol_grad = 1E-3
                    mol_opt_temp = optimize(mf, **params_pbe)
                    mf_opt_temp = UKS(mol_opt_temp, nxc=functionals[functional], nxc_kind='grid').density_fit()
                    mf_opt_temp.conv_tol = 1E-6
                    mf_opt_temp.conv_tol_grad = 1E-3
                    # mf_opt_temp.damp = 0.5  # if you need damping for better SCF convergence
                    # mf_opt_temp.max_cycle = 200 # increase SCF cycles if needed
                    mol_opt = optimize(mf_opt_temp, **conv_params_rel)
                    mf_opt = UKS(mol_opt, nxc=functionals[functional], nxc_kind='grid').density_fit()
                    # mf_opt.damp = 0.5 # uncomment if you reach convergence problems
                    mf_opt.max_cycle = 250
                    mf_opt.kernel()
                    -------------------------
                    """
                elif functionals[functional] == 'r2scan':
                    mf.xc = 'pbe'
                    mf.conv_tol = 1E-6
                    mf.conv_tol_grad = 1E-3
                    mol_opt_temp = optimize(mf, **params_pbe)
                    mf_opt_temp = dft.UKS(mol_opt_temp, xc=functionals[functional]).density_fit()
                    mf_opt_temp.conv_tol = 1E-6
                    mf_opt_temp.conv_tol_grad = 1E-3
                    # mf_opt_temp.damp = 0.5  # if you need damping for better SCF convergence
                    # mf_opt_temp.max_cycle = 250 # increase SCF cycles if needed
                    mol_opt = optimize(mf_opt_temp, **conv_params_rel)
                    mf_opt = dft.UKS(mol_opt, xc=functionals[functional]).density_fit()
                    # mf_opt.damp = 0.5 # uncomment if you reach convergence problems
                    mf_opt.max_cycle = 250
                    mf_opt.kernel()
                elif functionals[functional] == 'pbe0':
                    mol_opt = mol_pbe0_opt
                    mf_opt = dft.UKS(mol_opt, xc=functionals[functional]).density_fit()
                    # mf_opt.damp = 0.5 # uncomment if you reach convergence problems
                    mf_opt.max_cycle = 250
                    mf_opt.kernel()
                else:
                    mf.xc = functionals[functional]
                    mol_opt = optimize(mf, **conv_params)
                    mf_opt = dft.UKS(mol_opt, xc=functionals[functional]).density_fit()
                    # mf_opt.damp = 0.5 # uncomment if you reach convergence problems
                    mf_opt.max_cycle = 250
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
                if functionals[functional] in ['GGA_XC_PBE', 'MGGA_XC_SCAN', 'MGGA_HM', 'GGA_HM']:

                    """
                    !!!
                    Run this code with only one of these four functionals together
                    with the other functionals.
                    !!!
                    """

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
                    file.write("alpha-HOMO energy: {:.15f} Hartree\n".format(energy_alpha_homo))
                    file.write("alpha-LUMO energy: {:.15f} Hartree\n".format(energy_alpha_lumo))
                    file.write("beta-HOMO energy: {:.15f} Hartree\n".format(energy_beta_homo))
                    file.write("beta-LUMO energy: {:.15f} Hartree\n".format(energy_beta_lumo))
                    file.write("alpha-HOMO-LUMO gap: {:.15f} Hartree\n".format(alpha_hl_gap))
                    file.write("beta-HOMO-LUMO gap: {:.15f} Hartree\n".format(beta_hl_gap))


            # saving cube files for alpha/beta homo and lumo orbitals
            if not config.getboolean('system', 'single_atm'):
                save_cubefile_alpha_beta(mol_opt, mf_opt, output_folder, alpha_homo_idx, alpha_lumo_idx, beta_homo_idx, beta_lumo_idx)
                save_orbitals(mf_opt, output_folder, alpha_homo_idx, alpha_lumo_idx, beta_homo_idx, beta_lumo_idx)