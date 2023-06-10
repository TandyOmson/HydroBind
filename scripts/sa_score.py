#!/bin/python3

# Intermediate analysis script for SA score of hydrocarbons in results/Dhydrocarbons

import pandas as pd
import subprocess as sp
from rdkit import Chem
from rdkit.Chem import Draw

df = pd.read_pickle('/home/spine/DProjects/DhydroEn/results/Dhydrocarbons/all_mols_2.pkl')

df = df[df["mol"] != "InvalidSMILES"]

print(df)
df = df[df['sa_score'] > 3.0]
print(df['sa_score'].describe())

# Find molecule IDs with SA score < mean - 3*sdev
df = (df[df['sa_score'] < df['sa_score'].mean() - 3*df['sa_score'].std()])
df.sort_values(by=['sa_score'], inplace=True, ascending=True)
print(df)

def show_mol(mol):
    if isinstance(mol, Chem.rdchem.Mol):
        # 2D depiction
        #Draw.MolToFile(mol, "mol.png")
        #sp.run(["xdg-open", "mol.png"])
        # 3D depiction
        Chem.MolToPDBFile(mol, "mol.pdb")
        sp.Popen(["jmol","--silent","--nosplash","mol.pdb","&"])
        getc = input("Press ENTER to close: ")
        if getc == "":
            sp.run(["killall","java"])
        return 
    else:
        print("Invalid mol")
        return "InvalidSMILES"
    
getc = input("Show molecule? (y/n): ")
if getc == "y":
    for i in df.index:
        print(df.loc[i])
        show_mol(df.loc[i]["mol"])
        getc = input(f"Show next molecule? (y/n): ")
        if getc == "n":
            break