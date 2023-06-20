# Finds whether a mol object is an exclusion complex
# Calculates the centroid difference and number of atoms in the cavity
# Calls from miniball_container.cpp
# Contains the get_smallest_enclosing_sphere which calculates the radius around the centroid of the smallest enclosing sphere

# Import c++ extension container with pybind11
from miniball_container import get_smallest_enclosing_sphere
from mol_ops import xyz_to_mol
import numpy as np
from rdkit import Chem

def exclusion_check(complexmol,molId,centroidthreshold=4):
    """ Check whether the given mol object is an exclusion complex
    Exo complex if the difference between the centroid of the guest and host is > 4 A
    or the number of atoms in the cavity is < 6
    """
    # Get host and guest mol objects by splitting the complex (complex has already has been divided into residues)
    host_coords = []
    guest_coords = []

    for atom in complexmol.GetAtoms():
        if atom.GetPDBResidueInfo().GetResidueName() == "HOS":
            host_coords.append(complexmol.GetConformer().GetAtomPosition(atom.GetIdx()))
        elif atom.GetPDBResidueInfo().GetResidueName() == "GUE":
            guest_coords.append(complexmol.GetConformer().GetAtomPosition(atom.GetIdx()))

    # Miniball data consists of the centroid and radius squared of the smallest enclosing sphere
    # Radius squared is [-1], centroid is [0:3]
    miniball_data = get_smallest_enclosing_sphere(host_coords)

    # This gives smallest sphere which encloses the cavity, where a guest can fit
    host_radius = np.sqrt(miniball_data[-1])
    host_centroid = miniball_data[0:3]

    # Get guest centroid
    guest_coords = np.array(mol.GetConformer().GetPositions())
    guest_centroid = np.array(guest_coords).mean(axis=0)

    # Calculate distance between guest and host centroid
    centroid_diff = np.linalg.norm(guest_centroid - host_centroid)

    # Calculate number of atoms in cavity
    cavity_atoms = 0
    for atm in guest_coords:
        if np.linalg.norm(atm - host_centroid) < host_radius:
            cavity_atoms += 1

    # Check if exclusion complex
    isExo = False
    if centroid_diff > centroidthreshold or cavity_atoms < 6:
        isExo = True

    return isExo, centroid_diff, cavity_atoms

import sys

if __name__ == "__main__":
    mol = xyz_to_mol("example_complex.xyz",hostresname="HOS",guestresname="GUE")
    molId = 1
    print(exclusion_check(mol,molId))
