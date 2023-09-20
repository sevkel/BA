import numpy as np
from pyscf import gto, scf, cc


mol = gto.M(
    atom='H 0 0 0',  # in Angstrom
    basis='def2-qzvp',
    spin=1
)
mf = scf.HF(mol).run()
mycc = cc.CCSD(mf).run()
H2_en = mycc.e_tot


distances = np.arange(0.2,8,0.1)
energies=[]
for dist in distances:
    mol = gto.M(
        atom=f'H 0 0 0; H 0 0 {dist}',  # in Angstrom
        basis='def2-qzvp',
        symmetry=True,
    )
    mf = scf.HF(mol).run()
    mycc = cc.CCSD(mf).run()
    et = mycc.ccsd_t()
    energies.append(mycc.e_tot* 627.5096080305927-2*H2_en* 627.5096080305927)


destination = './resultsDFT/energies_H2_diss-longerrange'
np.savetxt(destination + '/CCSD.txt', energies)