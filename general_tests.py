import pyscf
from pyscf import dft, gto
# Definiere das Molekül
mol = pyscf.gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)

# Führe eine DFT-Rechnung durch
mf = pyscf.dft.RKS(mol)
mf.xc = 'lda,vwn'  # Hier kannst du den gewünschten Austausch-Korrelationsfunktional auswählen
mf.kernel()

# Erhalte die berechneten Molekülorbitale
mo_coeff = mf.mo_coeff
mo_energy = mf.mo_energy

# Bestimme die Anzahl der Elektronen im Molekül
num_electrons = mol.nelectron

# Bestimme die Energien der HOMO- und LUMO-Orbitale
homo_idx = num_electrons // 2 - 1
lumo_idx = homo_idx + 1
homo_energy = mo_energy[homo_idx]
lumo_energy = mo_energy[lumo_idx]

# Berechne die Homo-Lumo-Lücke
homo_lumo_gap = lumo_energy - homo_energy

# Gib die Ergebnisse aus
print("HOMO energy: {:.6f} Hartree".format(homo_energy))
print("LUMO energy: {:.6f} Hartree".format(lumo_energy))
print("HOMO-LUMO gap: {:.6f} Hartree".format(homo_lumo_gap))



