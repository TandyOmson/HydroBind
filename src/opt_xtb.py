""" Script for a generic optimisation using xTB, choosing whether to include solvent etc.
"""

import subprocess as sp
from rdkit import Chem

def opt(mol,molId,outfile,optfile,inp):
    """ optimises a given molecule using xTB """
    infile = f"{molId}.pdb"
    Chem.MolToPDBFile(mol,f"{infile}")
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

    guestmol = Chem.MolFromPDBFile(f"{optfile}",removeHs=False,sanitize=False)

    return guestmol
