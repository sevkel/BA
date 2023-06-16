import pyscf
from pyscf import dft, gto
import os
# Definiere das Molekül
mol = pyscf.gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'


mf = pyscf.dft.RKS(mol)
mf.xc = 'PBE'

e_tot = mf.kernel()
### Option 1: get_veff() returns J + V_xc
veff = mf.get_veff()
xc_energyone = -veff.exc
###

print("Austausch-Korrelationsenergie:", xc_energyone)


### Option 2: w and w/ XC part -> get the difference between the energies
mf.xc = ''
e_no_xc = mf.kernel()

# Berechne die Austausch-Korrelationsenergie
xc_energytwo = e_tot - e_no_xc

# Drucke die Austausch-Korrelationsenergie
print("Austausch-Korrelationsenergie:", xc_energytwo)
###


