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

def mol_gen(smiles):
    """in: smiles string
    - checks the smiles string for validity
    - returns a mol object to append to dictionary
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
            return mol
        except ValueError:
            return "InvalidSMILES"

def mol_out(mol, molId, **outdir):
    """in: molecule identifier (molid) and rdkit mol object
    - runs reduce to massage the structure and verify
    - returns a dictionary with the molid and the reduced mol object
    """
    print(f"Running reduce on {molId}")
    guestfile = str(molId) + "_rdkit" ".pdb"
    guestfile_reduced = str(molId) + "_reduced" + ".pdb"

    if outdir != None:
        guest_path = f"{outdir['outdir']}/{guestfile}"
        Chem.MolToPDBFile(mol, guest_path)
        reduced_path = f"{outdir['outdir']}/{guestfile_reduced}"
        sp.run(["reduce","-build",f"{outdir['outdir']}/{guestfile}"],stdout=open(reduced_path,"w"))
        mol = Chem.MolFromPDBFile(reduced_path,removeHs=False,sanitize=False)

        return mol

    else:
        sp.run(["reduce",f"{guestfile}"],stdout=open(f"{guestfile_reduced}","w"))
        mol = Chem.MolFromPDBFile(f"{guestfile_reduced}",removeHs=False,sanitize=False)
        sp.run(["rm",f"{guestfile}"])
    
        return mol

def single_smiles(smiles, molId, **outdir):
    """ Takes in a smiles string
    reutrns a dictionary with a reduced mol object
    """
    moldict = {'molId':molId}
    moldict['SMILES'] = smiles
    moldict['mol'] = mol_gen(moldict['SMILES'])
    moldict['path'] = None

    if moldict['mol'] != "InvalidSMILES":
        if outdir:
            moldict['mol'] = mol_out(moldict['mol'], moldict['molId'], outdir=outdir['outdir'])
        else:
            moldict['mol'] = mol_out(moldict['mol'], moldict['molId'])

    return moldict['mol']

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
