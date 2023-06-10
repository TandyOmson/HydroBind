# Uses xTB to dock generated molecules into the host molecule (CB7)

from rdkit import Chem
import subprocess as sp
import pandas as pd

from mol_ops import xyz_to_mol

def dock(mol,molId,hosttopo,xtbinp,**outdir):
    """ in: mol object, molId, host topology file, xTB input file, could add options for docking (ie. --pocket, --nostack, --noangular)
    Docks the mol object into the host topology
    returns a dictionary with the molId and the docked mol object
    optionally returns the xtb output and optimised topologies from guest and docking
    """
    
    # Optimize guest topology (same input file can be used for both optimization and docking)
    guesttopo = f"{molId}.pdb"
    Chem.MolToPDBFile(mol,guesttopo)
    if outdir:
        sp.run(["xtb","--input",f"{xtbinp}",f"{guesttopo}","--opt","--alpb","water"],stdout=open(f"{outdir['outdir']}/{molId}_opt.out","w"))
    else:
        sp.run(["xtb","--input",f"{xtbinp}",f"{guesttopo}","--opt","--alpb","water"])

    # Dock guest into host
    if outdir:
        sp.run(["xtb","dock",f"{hosttopo}","xtbopt.pdb","--input",f"{xtbinp}","--onlypocket"],stdout=open(f"{outdir['outdir']}/{molId}_dock.out","w"))
    else:
        sp.run(["xtb","dock",f"{hosttopo}","xtbopt.pdb","--input",f"{xtbinp}","--onlypocket"])
    
    # Read in docked guest
    complextopo = "pocket.xyz"
    mol = xyz_to_mol(complextopo)

    if outdir:
        sp.run(["cp","xtbopt.pdb",f"{outdir['outdir']}/{molId}_xtbopt.pdb"])
        Chem.MolToPDBFile(mol,f"{outdir['outdir']}/{molId}_docked.pdb")
    
    # Clean up
    sp.run(["rm","xtbopt.pdb"])
    sp.run(["rm","pocket.xyz"])
    sp.run(["rm",f"{molId}.pdb"])
    sp.run(["rm","xtbtopo.mol"])
    sp.run(["rm","wbo"])
    sp.run(["rm","xtbrestart"])
    sp.run(["rm","xtbopt.log"])
    sp.run(["rm",".xtboptok"])

    return mol

# To move to test
if __name__ == "__main__":
    # Grab imformation from DataFrame
    df = pd.read_pickle("/home/spine/DProjects/DhydroEn/results/Dhydrocarbons/all_mols_2.pkl")
    # Drop molecules that have a SA score < 3
    df = df[df["sa_score"] > 3]

    # xTB input file
    xtbinp = "/home/spine/DProjects/DhydroEn/data/xtbinp.inp"
    # Host topology file
    hosttopo = "/home/spine/DProjects/DhydroEn/data/DCB7/xtbopt.pdb"

    # Dock molecules
    mol = df["mol"].values
    molId = df.index.values

    z=0
    moldict = {}
    for i,j in zip(mol,molId):
        while z < 1:
            complexmol = dock(i,j,hosttopo,xtbinp,outdir="/home/spine/DProjects/DhydroEn/results/DdockedHydrocarbons")
            moldict[j] = complexmol
            z+=1
        
    complexdf = pd.DataFrame.from_dict(moldict,orient="index",columns=["docked_mol"])
    df = pd.concat([df,complexdf],axis=1)
    df.to_pickle("/home/spine/DProjects/DhydroEn/results/Dhydrocarbons/all_mols_3.pkl")
