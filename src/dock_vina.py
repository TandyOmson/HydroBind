""" Module for docking using AutoDock Vina
Inputs: Molecule object, Molecule ID, Host topology file, Vina input file
Outputs: Docking output file, and file with best pose
Returns: Complex molecule object
"""

import subprocess as sp
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem

# Import tblite-python module, a python API for singlepoint tight binding calculations (a lightweight energy evaluation)
#from tblite.interface import Calculator
#import numpy as np

from mol_ops import change_complex_resnames, ammend_pdb_spacing

def dock(mol,molId,hostfile,inp):
    """ Dock a molecule into a host
    hostfile is a .pdbqt file, use command "obabel host.pdb -O host.pdbqt -xrh" to convert
    """
    
    infile = f"{molId}.pdb"
    pdbqtfile = f"{molId}.pdbqt"

    dockoutfile = f"{molId}_poses.pdb"
    posepdbqtfile = f"{molId}_poses.pdbqt"
    posefile = f"{molId}_complex.pdb"

    Chem.MolToPDBFile(mol,f"{infile}")

    # Convert guest .pdb to .pdbqt
    sp.run(["obabel",f"{infile}","-O",f"{pdbqtfile}","-xh"],stderr=open("obabel.log","w"))

    # Run Vina
    sp.run(["vina",f"--config={inp}",f"--receptor={hostfile}",f"--ligand={pdbqtfile}",f"--out={posepdbqtfile}"],stdout=open("vina.log","w"))

    # Revert output .pdbqt to .pdb
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

    # TESTINGPHASE: Brief optimisation of all conformers
    #complexmols.UpdatePropertyCache(strict=False)
    #Chem.GetSymmSSSR(complexmols)
    #complexmols.GetRingInfo().NumRings()
    #
    #AllChem.MMFFOptimizeMoleculeConfs(complexmols, ignoreInterfragInteractions=False, nonBondedThresh=100.0)

    # TESTINGPHASE: Quick GNF2 singlepoint energy evaluation to score poses
    # Set as a property of each conformer in complexmols
    #for i in range(num_confs):
    #    mol = complexmols.GetConformer(i)
    #    atomicnums = np.array([atom.GetAtomicNum() for atom in complexmols.GetAtoms()])
    #    coords = np.array(mol.GetPositions())
    #
    #    calc = Calculator("GFN2-xTB", atomicnums, coords)
    #    res = calc.singlepoint()
    #    print(res.get("energy"))
    #
    #    complexmols.GetConformer(i).SetProp("energy",str(res.get("energy")))
    #
    #for c in complexmols.GetConformers():
    #    print(c.GetProp("energy"))

    # Take only the best pose
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
    sp.run(["rm","vina.log"])
    sp.run(["rm","obabel.log"])

    return complexmol

if __name__ == "__main__":
    guestmol = Chem.MolFromPDBFile('7368.pdb',removeHs=False,sanitize=False)
    inp = "vinadock.inp"
    complexmol = dock(guestmol,7368,"host.pdbqt",inp)
