""" Docking method that aligns the host and guest axis, places the coordinates of the guest inside the host in a few configurations, and optimises
    Takes the best of a few binding poses
    This method uses PCA (principal component analysis) to align the host and guest axis
    In short, PCA finds the direction of the largest variance in the data, and aligns the data to that direction
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.decomposition import PCA
import numpy as np
from mol_ops import xyz_to_mol

def dock(mol,molID,hosttopo,dockoutfile,posefile):

    # Align host and guest axis
    # Get coordinates of guest
    mol = Chem.rdmolops.RemoveHs(mol)

    guest_atoms = [atm.GetSymbol() for atm in mol.GetAtoms()]
    guest_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for atm in mol.GetAtoms()])

    hostmol = Chem.MolFromPDBFile(hosttopo,removeHs=False,sanitize=False)
    host_coords = np.array([hostmol.GetConformer().GetAtomPosition(atm.GetIdx()) for atm in hostmol.GetAtoms() if atm.GetAtomicNum() != 1])

    #add clouds of points around the atoms of the guest molecule, 3 sections in polar (theta) and 6 in azimuthal (phi)
    #This prepares the coordinates for the PCA
    #   x = r*sin(theta)*cos(phi)
	#	y = r*sin(theta)*sin(phi)        
	#	z = r*cos(theta)
    atomic_radii = {'H' :1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'Cl':1.75, 'Br':1.85, 'I':1.98, 'P':1.8, 'S':1.8, 'As':1.85, 'B':2.13, 'Si':2.1, 'Se':1.9, 'Te':2.06}
    
    cloud_points = [] 
    for index,atm in enumerate(guest_atoms):
        cloud_radius = atomic_radii[atm]
        # Add a cloud point above and below the atom
        cloud_points.append(guest_coords[index] + np.array([0,0,cloud_radius]))
        cloud_points.append(guest_coords[index] + np.array([0,0,-cloud_radius]))
        # Add a cloud of points around the atom, 3 in theta and 6 in phi
        for theta in np.linspace(0,np.pi,3):
            for phi in np.linspace(0,2*np.pi,6):
                cloud_points.append(guest_coords[index] + np.array([cloud_radius*np.sin(theta)*np.cos(phi),cloud_radius*np.sin(theta)*np.sin(phi),cloud_radius*np.cos(theta)]))

    # Add cloud points to the guest coordinates
    cloud_points = np.array(cloud_points)
    guest_coords_with_clouds = np.concatenate((guest_coords,cloud_points),axis=0)

    # Initiliase PCA
    pca = PCA(n_components=3)
    # Fit PCA to guest coordinates and transform the coordinates
    pca.fit_transform(guest_coords_with_clouds)
    transform_coord = pca.transform(guest_coords)

    # Centre the transformed coordinates on the host centroid
    host_centroid = np.mean(host_coords,axis=0)
    transform_coord_centred = transform_coord - np.mean(transform_coord,axis=0) + host_centroid

    # Set molecule coordinates to the transformed coordinates
    for index,atm in enumerate(mol.GetAtoms()):
            mol.GetConformer().SetAtomPosition(atm.GetIdx(),transform_coord_centred[index])
            
    # Add Hs back to the guest
    mol = Chem.AddHs(mol)

    # Finally combine the host and guest coordinates
    complexmol = Chem.CombineMols(hostmol,mol)

    # Converge the molecule
    n_steps = 1000000
    tol=1e-8
    AllChem.MMFFSanitizeMolecule(complexmol)
    ff = AllChem.MMFFGetMoleculeForceField(complexmol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(complexmol, mmffVariant='MMFF94', mmffVerbosity = 0), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
    ff.Initialize()
    cf = ff.Minimize(maxIts=n_steps,energyTol=tol,forceTol=tol)
    en = ff.CalcEnergy()

    # Write out the docked molecule
    Chem.MolToXYZFile(complexmol,dockoutfile)

    # Read in docked guest from .xyz file with formatted residues
    mol = xyz_to_mol(dockoutfile)

    return mol

if __name__ == "__main__":
    guestmol = Chem.MolFromPDBFile('7368_opt.pdb',removeHs=False,sanitize=False)
    complex = dock(guestmol,1,"host.pdb","best.xyz","pose.xyz")
