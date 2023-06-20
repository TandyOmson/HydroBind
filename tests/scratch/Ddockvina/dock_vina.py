""" Module for docking using AutoDock Vina
Inputs: Molecule object, Molecule ID, Host topology file, Vina input file
Outputs: Docking output file, and file with best pose
Returns: Complex molecule object
"""

import subprocess as sp
from rdkit import Chem
from rdkit import RDLogger
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
    sp.run(["vina",f"--config={inp}",f"--ligand={pdbqtfile}",f"--out={posepdbqtfile}"],stdout=open("vina.log","w"))

    # Revert output .pdbqt to .pdb
    sp.run(["obabel",f"{posepdbqtfile}","-O",f"{dockoutfile}"],stderr=open("obabel.log","a"))

    # Disable RDKit logging
    RDLogger.DisableLog('rdApp.*')

    # Extract best conformer (make it only conformer to be compatible with Chem.CombineMols)
    guestmols = Chem.MolFromPDBFile(f"{dockoutfile}",removeHs=False,sanitize=False)
    guestmol = Chem.Mol(guestmols)
    guestmol.RemoveAllConformers()
    guestmol.AddConformer(guestmols.GetConformer(),assignId=True)

    # Combine host and guest
    hostmol = Chem.MolFromPDBFile(f"{hostfile}",removeHs=False,sanitize=False)
    complexmol = Chem.CombineMols(hostmol,guestmol)

    # OPTIONAL: Optimse complex (this could be implemented for all conformers as well
    # requires hostmol to have the same number of copies of conformers as conformers in guestmols before combinemols
    #complexmol.UpdatePropertyCache(strict=False)
    #Chem.GetSymmSSSR(complexmol)
    #complexmol.GetRingInfo().NumRings()

    #AllChem.MMFFOptimizeMolecule(complexmol, ignoreInterfragInteractions=False, nonBondedThresh=100.0)

    # Write out docked molecule
    complexmol = change_complex_resnames(complexmol,"GUE","HOS")
    Chem.MolToPDBFile(complexmol,posefile)
    ammend_pdb_spacing(posefile)

    # Clean up files not standard to class method
    sp.run(["mv",f"{dockoutfile}",f"{molId}_dock.out"])
    sp.run(["rm",f"{pdbqtfile}"])
    sp.run(["rm",f"{posepdbqtfile}"])
    sp.run(["rm","vina.log"])
    sp.run(["rm","obabel.log"])

    return complexmol

if __name__ == "__main__":
    guestmol = Chem.MolFromPDBFile('7368.pdb',removeHs=False,sanitize=False)
    inp = "vinadock.inp"
    complexmol = dock(guestmol,7368,"host.pdbqt",inp)