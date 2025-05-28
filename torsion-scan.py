from ase.io import read, write
from ase.atoms import Atoms
from ase.neighborlist import neighbor_list
import numpy as np
import sys
import os
from collections import deque

def find_rotatable_atoms(mol, fixed_atom, start_atom):
    """BFS to find all atoms connected to start_atom, avoiding traversal back to fixed_atom."""
    graph = neighbor_list('i', mol, cutoff=1.6), neighbor_list('j', mol, cutoff=1.6)
    bonds = [[] for _ in range(len(mol))]
    for i, j in zip(*graph):
        bonds[i].append(j)
        bonds[j].append(i)

    visited = set()
    queue = deque([start_atom])
    visited.add(fixed_atom) 
    rotatable = set()

    while queue:
        current = queue.popleft()
        rotatable.add(current)
        for neighbor in bonds[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    return rotatable

# Atom indices are 0-based (Python convention)
xyz_file = sys.argv[1]
i, j, k, l = map(int, sys.argv[2:6])
mol = read(xyz_file)
# Define scan angles (in degrees)
angles = np.linspace(0, 180, 37, endpoint=True)
base_name = os.path.splitext(os.path.basename(xyz_file))[0]
# Get which atoms should rotate
rotatable_atoms = find_rotatable_atoms(mol, k, l)
mask = [1 if idx in rotatable_atoms else 0 for idx in range(len(mol))]
for idx, angle in enumerate(angles):
    mol.set_dihedral(i,j,k,l,angle,mask=mask)
    write(f"{base_name}_scan_{int(angle)}.xyz", mol)
print(f"Generated {len(angles)} torsion conformers.")



