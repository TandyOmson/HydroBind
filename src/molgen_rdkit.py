# conformer generator - see a short manual version in scripts
# ----------------------------
# Intake SMILES string, each molecule has a pd.Series
# Use RDkit to generate a conformer
# Use AMBER's reduce to massage the structure and verfiy
# Rapid conformer search to get a reasonable guess at final geometry
# ----------------------------

# Depedenices
# ----------------------------
# RDkit
# AmberTools
# ----------------------------

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import subprocess as sp

def mol_gen(smiles, molId):
    """in: smiles string
    - checks the smiles string for validity
    - returns a mol object to append to dictionary
    - runs reduce to massage the structure and verify presence of all hydrogens
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "InvalidSMILES"
    else:
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol,useRandomCoords=True,randomSeed=0xf00d,ignoreSmoothingFailures=True,randNegEig=True)
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except ValueError:
            return "InvalidSMILES"

    guestfile = f"{molId}_rdkit.pdb"
    guestfile_reduced = f"{molId}_reduced.pdb"

    Chem.MolToPDBFile(mol, guestfile)
    sp.run(["reduce","-build",f"{guestfile}"],stdout=open(guestfile_reduced,"w"))
    mol = Chem.MolFromPDBFile(guestfile_reduced,removeHs=False,sanitize=False)

    return mol

# To move to test
if __name__ == "__main__":
    smipath = "./../data/hydroIDs.smi"
    smi_gen = (i.split()[1] for i in open(smipath,"r").readlines()[1:])
    all_mols = pd.DataFrame()
    
    for i in range(8992):
        moldict = single_smiles(next(smi_gen),i,outdir="./../results/Dhydrocarbons")
        all_mols = pd.concat([all_mols,pd.DataFrame([moldict])],ignore_index=True)

    all_mols = all_mols.set_index("molId")
    print(all_mols)
    all_mols.to_pickle("./../results/Dhydrocarbons/all_mols.pkl")
