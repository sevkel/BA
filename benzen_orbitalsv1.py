from pyscf import gto, tools, dft, scf
import numpy as np

mol = gto.Mole()
mol.atom = 'C6H6.xyz'
mol.basis = 'sto-3g'
mol.build()

mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()

# Extrahieren der Orbitale
energies = mf.mo_energy
orbitals = mf.mo_coeff

# Finde das HOMO-Orbital; ENERGIES ARE ALREADY SORTED!
#homo_index = np.where(energies <= 0)[0][-1]
#homo_index = 21
#homo_orbital = orbitals[homo_index]

# Finde das LUMO-Orbital
lumo_idx = mf.mo_occ.tolist().index(0.0)
lumo_orbital = orbitals[lumo_idx]

homo_idx = lumo_idx - 1
homo_orbital = orbitals[homo_idx]
dateiname= "orbitals.txt"

# Erstellen der Cube-Datei für das homo Orbital
#tools.cubegen.orbital(mol, 'orbital_homo3.cube', orbitals[:, homo_idx])

# Erstellen der Cube-Datei für das lumo Orbital
#tools.cubegen.orbital(mol, "orbital_lumo3.cube", orbitals[:, lumo_idx])


#Ausgabe der Ergebnisse
#print("energies:\n", energies)
#print("HOMO-energy:\n", energies[homo_idx])
'''
print("Orbitale\n", orbitals)

print("HOMO-Orbital:")
print(homo_orbital)

print("\nLUMO-Orbital:")
print(lumo_orbital)
'''
#print("Homo-Index:", homo_idx)
#print("Lumo-Index:", lumo_idx)

#mf.analyze() #IMPORTANT
print(mol.atom_coords(unit='ang'))

from pyscf.geomopt.geometric_solver import optimize