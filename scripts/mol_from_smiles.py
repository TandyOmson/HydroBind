# USAGE: python mol_from_smiles.py <SMILES_STRING> <OUTPUT_FILE>
# Example: python mol_from_smiles.py "CC(C)(C)C1=CC=C(C=C1)C(=O)O" "test.pdb"

# Conformer generator
# ----------------------------
# Intake SMILES string, each molecule has a pd.Series
# Use RDkit to generate a conformer
# Use AMBER's reduce to massage the structure and verfiy
# ----------------------------

# Depedenices
# ----------------------------
# rdkit
# AmberTools23
# ----------------------------

from rdkit import Chem
from rdkit.Chem import AllChem
import sys

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

def mol_out(mol,outfile):
    """in: molecule identifier (molid) and rdkit mol object
    - runs reduce to massage the structure and verify
    - returns a dictionary with the molid and the reduced mol object
    """
    if mol != "InvalidSMILES":
        Chem.MolToPDBFile(mol, outfile)
    else:
        print("Not a valid SMILES string")

    return

# Simple script usage
smile = sys.argv[1]
molfile= sys.argv[2]

mol_out(mol_gen(smile),molfile)