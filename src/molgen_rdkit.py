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
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign

import pandas as pd

def mol_gen(smiles, molId):
    """in: smiles string
    - checks the smiles string for validity
    - returns a mol object to append to dictionary
    - runs reduce to massage the structure and verify presence of all hydrogens
    """

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))

    # Check for validity
    param = rdDistGeom.ETKDGv2()
    param.pruneRmsThresh = 0.2
    rdDistGeom.EmbedMolecule(mol,params=param)

    AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s")

    # Write to pdb
    Chem.MolToPDBFile(mol,f"{molId}_rdkit.pdb")

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
