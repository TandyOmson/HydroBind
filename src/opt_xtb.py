""" Script for a generic optimisation using xTB, choosing whether to include solvent etc.
"""

import subprocess as sp
from rdkit import Chem
from mol_ops import ammend_pdb_spacing

def opt(mol,molId,inp):
    """ optimises a given molecule using xTB """
    infile = f"{molId}.pdb"
    outfile = f"{molId}_opt.out"
    optfile = f"{molId}_opt.pdb"

    Chem.MolToPDBFile(mol,f"{infile}")
    # MolToPDBFile formats differently for complexes and single molecules, this is a temporary fix
    if len(mol.GetAtoms()) > 125:
        ammend_pdb_spacing(f"{infile}")

    print("Command: xtb --input",f"{inp}",f"{infile}","--opt","--alpb","water",">",f"{outfile}")
    sp.run(["xtb","--input",f"{inp}",f"{infile}","--opt","--alpb","water"],stdout=open(outfile,"w"))
    sp.run(["mv","xtbopt.pdb",f"{optfile}"])

    # Clean up
    sp.run(["rm",f"{infile}"])
    sp.run(["rm","charges"])
    sp.run(["rm","xtbtopo.mol"])
    sp.run(["rm","wbo"])
    sp.run(["rm","xtbrestart"])
    sp.run(["rm","xtbopt.log"])
    sp.run(["rm",".xtboptok"])

    mol = Chem.MolFromPDBFile(f"{optfile}",removeHs=False,sanitize=False)

    return mol
