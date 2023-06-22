""" Docking method for guests in CB7 
Three step method:
- Step 1: AutoDock Vina using large search box
- Step 2: AutoDock Vina using smaller search box to fit inner cavity
- Step 3: Use alignment docking in an edge case
"""

from vina import Vina
from rdkit import Chem
import subprocess as sp
from rdkit import RDLogger
import numpy as np

from mol_ops import change_complex_resnames, ammend_pdb_spacing

def dock(mol,molId,hostfile,inp):
    """ Dock a molecule into a host
    hostfile is a .pdbqt file, use command "obabel host.pdb -O host.pdbqt -xrh" to convert
    """
    
    infile = f"{molId}.pdb"
    pdbqtfile = f"{molId}.pdbqt"

    posepdbqtfile = f"{molId}_poses.pdbqt"
    dockoutfile = f"{molId}_poses.pdb"
    posefile = f"{molId}_complex.pdb"

    Chem.MolToPDBFile(mol,f"{infile}")

    # Convert guest .pdb to .pdbqt
    sp.run(["obabel",f"{infile}","-O",f"{pdbqtfile}","-xh"],stderr=open("obabel.log","w"))

    # Run Vina with large box
    v = Vina(sf_name='vina')
    v.set_receptor(rigid_pdbqt_filename=hostfile)
    v.set_ligand_from_file(pdbqtfile)
    
    v.compute_vina_maps(center=[0,0,0], box_size=[24,24,24])
    v.optimize()

    # Exhaustiveness deterimes the number of MC samples
    v.dock(exhaustiveness=20, n_poses=10)
    big_box_ens = v.energies(n_poses=10)
    v.write_poses(posepdbqtfile[:-6]+"_big.pdbqt",overwrite=True)

    # Run Vina with small box
    v.compute_vina_maps(center=[0,0,0], box_size=[10.780,10.051,10.131])
    v.optimize()
    v.dock(exhaustiveness=20, n_poses=10)
    small_box_ens = v.energies(n_poses=10)
    v.write_poses(posepdbqtfile[:-6]+"_small.pdbqt",overwrite=True)

    # Revert output .pdbqt to .pdb
    sp.run(["cat",f"{posepdbqtfile[:-6]}_big.pdbqt",f"{posepdbqtfile[:-6]}_small.pdbqt"],stdout=open(f"{posepdbqtfile}","w"))
    sp.run(["obabel",f"{posepdbqtfile}","-O",f"{dockoutfile}"],stderr=open("obabel.log","a"))

    # Disable RDKit logging
    RDLogger.DisableLog('rdApp.*')

    # Read in all binding poses
    guestmols = Chem.MolFromPDBFile(f"{dockoutfile}",removeHs=False,sanitize=False)
    hostmol = Chem.MolFromPDBFile(f"{hostfile}",removeHs=False,sanitize=False)

    # Rewrite dockoutfile with hosts included
    num_confs = guestmols.GetNumConformers()
    # Make num_confs copies of the conformers in hostmol
    hostmols = Chem.Mol(hostmol)
    for i in range(num_confs-1):
        hostmols.AddConformer(hostmol.GetConformer(),assignId=True)

    # Combine host and guest
    complexmols = Chem.CombineMols(hostmols,guestmols)

    # Write out docked molecules in all binding poses
    Chem.MolToPDBFile(complexmols,dockoutfile)

    # Take the best pose from the large or small box
    # If the smaller box retruns an energy within 3 kcal/mol of the larger box, use the smaller box
    if small_box_ens[:,0].min() < big_box_ens[:,0].min() + 3:
        sp.run(["obabel",f"{posepdbqtfile[:-6]}_small.pdbqt","-O",f"{posefile}"],stderr=open("obabel.log","a"))
        guestmols = Chem.Mol(Chem.MolFromPDBFile(f"{posefile}",removeHs=False,sanitize=False))
        
    guestmol = Chem.Mol(guestmols)
    guestmol.RemoveAllConformers()
    guestmol.AddConformer(guestmols.GetConformer(),assignId=True)

    # Combine host and guest
    complexmol = Chem.CombineMols(hostmol,guestmol)

    # Write out docked molecule
    complexmol = change_complex_resnames(complexmol,"GUE","HOS")
    Chem.MolToPDBFile(complexmol,posefile)
    ammend_pdb_spacing(posefile)

    # Clean up files not standard to class method
    sp.run(["mv",f"{dockoutfile}",f"{molId}_dock.out"])
    sp.run(["rm",f"{infile}"])
    sp.run(["rm",f"{pdbqtfile}"])
    sp.run(["rm",f"{posepdbqtfile}"])
    sp.run(["rm",f"{posepdbqtfile[:-6]}_big.pdbqt"])
    sp.run(["rm",f"{posepdbqtfile[:-6]}_small.pdbqt"])
    sp.run(["rm","obabel.log"])

    return complexmol

if __name__ == "__main__":
    import time
    start = time.time()
    hostfile = "/home/spine/DProjects/DhydroEn/data/DCB7/host.pdbqt"
    guestmol = Chem.MolFromPDBFile("7368.pdb",removeHs=False,sanitize=False)
    inp = None
    complexmol = dock(guestmol,7368,hostfile,inp)
    end = time.time()
    print("time: ",end-start,"s")
