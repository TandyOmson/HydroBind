# Finds whether a mol object is an exclusion complex
# Defines exclusion complexes as either:
# Distance between centroid of host and guest > centroid_diff_threshold
# Number of atoms in the cavity < cavity_atoms_threshold

from mol_ops import xyz_to_mol
import numpy as np
from scipy.spatial import Delaunay
from rdkit import Chem

def exclusion_check(complexmol,molId,centroiddiffthreshold=4,cavityatomsthreshold=6):
    """ Check whether the given mol object is an exclusion complex
    Exo complex if the difference between the centroid of the guest and host is > 4 A
    or the number of atoms in the cavity is < 6
    """
    # Get host and guest coordinates by splitting the complex that has been formatted by xyz_to_mol
    host_coords = []
    guest_coords = []

    # Retrive host and guest coordinates, ignoring hydrogens
    for atom in complexmol.GetAtoms():
        if atom.GetPDBResidueInfo().GetResidueName() == "HOS":
            if atom.GetAtomicNum() != 1:
                host_coords.append(complexmol.GetConformer().GetAtomPosition(atom.GetIdx()))
        elif atom.GetPDBResidueInfo().GetResidueName() == "GUE":
            if atom.GetAtomicNum() != 1:
                guest_coords.append(complexmol.GetConformer().GetAtomPosition(atom.GetIdx()))

    # Get host and guest centroid
    guest_centroid = np.array(guest_coords).mean(axis=0)
    host_centroid = np.array(host_coords).mean(axis=0)

    # Calculate distance between guest and host centroid
    centroid_diff = np.linalg.norm(guest_centroid - host_centroid)

    # Delauny defines the convex hull of the host atoms
    # Delauny is a triangulation such that none of the host atoms are inside the circumsphere of any tetrahedron in the traingulation
    hull = Delaunay(host_coords)

    # Calculate number of atoms in cavity
    cavity_atoms = 0
    for atm in guest_coords:
        # Points outside the triangulation return -1
        if hull.find_simplex(atm) >= 0:
            cavity_atoms += 1

    # Check if exclusion complex
    isExo = False
    if centroid_diff > centroiddiffthreshold or cavity_atoms < cavityatomsthreshold:
        isExo = True

    return isExo, centroid_diff, cavity_atoms

import sys

if __name__ == "__main__":
    mol = xyz_to_mol("example_complex.xyz",hostresname="HOS",guestresname="GUE")
    molId = 1
    print(exclusion_check(mol,molId))
