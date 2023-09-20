import numpy as np
from pyscf import gto
from pyscf import dft
import density_functional_approximation_dm21 as dm21
from pylibnxc.pyscf import UKS


def compute_dissociation(bond_lengths, en_h_atom, func):
    energies = []

    for bond_length in bond_lengths:

        mol = gto.Mole()
        mol.atom = f'H 0 0 0; H 0 0 {bond_length}'
        mol.basis = 'def2-qzvp'
        mol.symmetry = 'Dooh'
        #mol.verbose = 4
        mol.build()
        if (func in ['DM21', 'DM21m', 'DM21mu', 'DM21mc']):
            mf = dft.UKS(mol)
            mf.xc = 'B3LYP'
            mf.run()
            dm0 = mf.make_rdm1()
            irrep_nelec0 = mf.get_irrep_nelec()

            mf = dft.UKS(mol)
            if func == 'DM21':
                mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
            elif func == 'DM21m':
                mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21m)
            elif func == 'DM21mu':
                mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mu)
            elif func == 'DM21mc':
                mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mc)
            mf.grids.level = 5
            mf.conv_tol = 1E-6
            mf.conv_tol_grad = 1E-3
            mf.irrep_nelec = irrep_nelec0
            mf.verbose = 4

            mf.kernel(dm0=dm0)

        elif (func in ['GGA_XC_PBE', 'MGGA_XC_SCAN', 'GGA_HM', 'MGGA_HM']):
            mf = UKS(mol, nxc=func, nxc_kind='grid')
            mf.grids.level = 5
            #mf.damp = 0.5
            mf.conv_tol = 1E-6
            mf.conv_tol_grad = 1E-3
            mf.kernel()
        else:
            mf = dft.UKS(mol)
            mf.xc = func
            mf.grids.level = 5
            #mf.damp = 0.5
            mf.conv_tol = 1E-6
            mf.conv_tol_grad = 1E-3
            mf.kernel()

        energy_H2 = mf.e_tot

        energies.append(energy_H2 * 627.5096080305927 - 2 * en_h_atom * 627.5096080305927)

    return energies

def compute_single_atom_energy(func):
    mol = gto.Mole()
    mol.atom = 'H 0 0 0'
    mol.basis = 'def2-qzvp'
    mol.spin = 1
    mol.build()

    if (func in ['DM21', 'DM21m', 'DM21mu', 'DM21mc']):
        mf = dft.UKS(mol)
        if func == 'DM21':
            mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21)
        elif func == 'DM21m':
            mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21m)
        elif func == 'DM21mu':
            mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mu)
        elif func == 'DM21mc':
            mf._numint = dm21.NeuralNumInt(dm21.Functional.DM21mc)
        mf.conv_tol = 1E-6
        mf.conv_tol_grad = 1E-3
        mf.kernel()
    elif (func in ['GGA_XC_PBE', 'MGGA_XC_SCAN', 'GGA_HM', 'MGGA_HM']):
        mf = UKS(mol, nxc=func, nxc_kind='grid')
        mf.conv_tol = 1E-6
        mf.conv_tol_grad = 1E-3
        mf.kernel()
    else:
        mf = dft.UKS(mol)
        mf.xc = func
        mf.conv_tol = 1E-6
        mf.conv_tol_grad = 1E-3
        mf.kernel()

    return mf.e_tot

def main():
    functionals = ['MGGA_HM']  # set your functionals you want to test

    bond_lengths = np.arange(0.2, 8, 0.1)

    for func in functionals:

        h_en = compute_single_atom_energy(func)

        energies = compute_dissociation(bond_lengths, h_en, func)

        destination = './resultsDFT/H2_diss/energies_H2_diss-longerrange'
        np.savetxt(destination + f'/{func}.txt', energies)


if __name__ == '__main__':
    main()


