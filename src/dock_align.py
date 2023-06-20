""" Docking method that aligns the host and guest axis, places the coordinates of the guest inside the host in a few configurations, and optimises
    Takes the best of a few binding poses
    This method uses PCA (principal component analysis) to align the host and guest axis
    In short, PCA finds the direction of the largest variance in the data, and aligns the data to that direction
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from sklearn.decomposition import PCA
import numpy as np
from mol_ops import change_complex_resnames, ammend_pdb_spacing

from prog_ops import coords_check

def dock(mol,molId,hostfile,inp):
    """ Dock a molecule into a host"""

    dockoutfile = f"{molId}_dock.out"
    posefile = f"{molId}_complex.pdb"

    guest_atoms = [atm.GetSymbol() for atm in mol.GetAtoms()]
    guest_coords = np.array([mol.GetConformer().GetAtomPosition(atm.GetIdx()) for atm in mol.GetAtoms()])

    hostmol = Chem.MolFromPDBFile(hostfile,removeHs=False,sanitize=False)
    host_atoms = [atm.GetSymbol() for atm in hostmol.GetAtoms()]
    host_coords = np.array([hostmol.GetConformer().GetAtomPosition(atm.GetIdx()) for atm in hostmol.GetAtoms()])

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

    # Direct the principal axis of the guest molecule towards the z-axis (the axis pointing through the cavity of the host)
    theta = np.arctan2(transform_coord[0,0],transform_coord[0,2])
    rotation_matrix = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]])
    transform_coord = np.matmul(rotation_matrix,transform_coord.T).T

    # Centre the transformed coordinates on the host centroid
    transform_coord_centered = transform_coord.copy()
    transform_coord_centered[:,0] = transform_coord[:,0] - np.mean(transform_coord[:,0])
    transform_coord_centered[:,1] = transform_coord[:,1] - np.mean(transform_coord[:,1])
    transform_coord_centered[:,2] = transform_coord[:,2] - np.mean(transform_coord[:,2])
    
    # Set molecule coordinates to the transformed coordinates
    for index,atm in enumerate(mol.GetAtoms()):
            mol.GetConformer().SetAtomPosition(atm.GetIdx(),transform_coord_centered[index])

    # Thanks to the symmetric structure of CBs, the only way to change docking pose is a slight rotation around the z-axis if more conformers are required

    # Finally combine the host and guest coordinates
    complexmolrdkit = Chem.CombineMols(hostmol,mol)

    complexmolrdkit.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(complexmolrdkit)
    complexmolrdkit.GetRingInfo().NumRings()

    # Turn off addHs warnings
    RDLogger.DisableLog('rdApp.*')
    AllChem.MMFFOptimizeMolecule(complexmolrdkit, ignoreInterfragInteractions=False, nonBondedThresh=100.0)
    with open(dockoutfile,'w') as f:
         f.write("NO DOCKING INFO FOR ALIGN METHOD")

    #with open(posefile,'w') as f:
    #    f.write(f"{len(transform_coord_centered)+len(host_coords)}\n\n")
    #    for atom,specie in zip(host_coords,host_atoms):
    #        f.write(f"{specie} {atom[0]} {atom[1]} {atom[2]}\n")
    #    for atom,specie in zip(transform_coord_centered,guest_atoms):
    #        f.write(f"{specie} {atom[0]} {atom[1]} {atom[2]}\n")

            
    #complexmol = Chem.MolFromXYZFile(posefile)

    # Converge the molecule
    #n_steps = 1000000
    #tol=1e-8
    #AllChem.MMFFSanitizeMolecule(complexmol)
    #ff = AllChem.MMFFGetMoleculeForceField(complexmol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(complexmol, mmffVariant='MMFF94', mmffVerbosity = 0), ignoreInterfragInteractions=False, nonBondedThresh=100.0)
    #ff.Initialize()
    #cf = ff.Minimize(maxIts=n_steps,energyTol=tol,forceTol=tol)
    #en = ff.CalcEnergy()
    #with open(dockoutfile,'w') as f:
    #    f.write(f"{cf}\n")
    #    f.write(f"{en}\n")
    
    # Write out docked molecule
    complexmolrdkit = change_complex_resnames(complexmolrdkit,"GUE","HOS")
    Chem.MolToPDBFile(complexmolrdkit,posefile)
    ammend_pdb_spacing(posefile)

    return complexmolrdkit

if __name__ == "__main__":
    guestmol = Chem.MolFromPDBFile('7368_opt.pdb',removeHs=False,sanitize=False)
    complex = dock(guestmol,1,"host.pdb",inp)
