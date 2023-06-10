# Uses xTB to dock generated molecules into the host molecule (CB7)

from rdkit import Chem
import subprocess as sp
import pandas as pd

from mol_ops import xyz_to_mol

def dock(mol,molId,hosttopo,dockoutfile,posefile,inp):
    """ in: mol object, molId, host topology file, xTB input file, could add options for docking (ie. --pocket, --nostack, --noangular)
    Docks the mol object into the host topology
    returns complex mol object
    """
    
    # Write guest topology to file
    guesttopo = f"{molId}.pdb"
    Chem.MolToPDBFile(mol,guesttopo)
    
    # Dock guest into host
    sp.run(["xtb","dock",f"{hosttopo}",f"{guesttopo}","--input",f"{inp}","--pocket"],stdout=open(f"{dockoutfile}","w"))

    # Read in docked guest
    complextopo = "best.xyz"
    mol = xyz_to_mol(complextopo)

    sp.run(["mv",f"{complextopo}",f"{posefile}"])

    sp.run(["rm","pocket.xyz"])
    sp.run(["rm",f"{molId}.pdb"])

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
