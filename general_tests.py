import pyscf
from pyscf import dft, gto
# Definiere das Molekül
mol = pyscf.gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)

# Führe eine DFT-Rechnung durch
mf = pyscf.dft.RKS(mol)
mf.xc = 'pbe'  # Hier kannst du den gewünschten Austausch-Korrelationsfunktional auswählen
mf.kernel()





