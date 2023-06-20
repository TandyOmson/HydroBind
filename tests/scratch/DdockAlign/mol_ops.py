# Extra molecule operations for import into the modules

from rdkit import Chem

CB_ATOM_NUM = 126

def xyz_to_mol(xyzfile,**resnames):
    """ Goes from the xyz file of a complex to a pdb file with residue information
    """

    mol = Chem.MolFromXYZFile(xyzfile)

    if len(resnames)!=2:
        resnames = {'hostresname':'UNL','guestresname':'UNL'}
        print("Residue names not specified. Using default names: UNL for both host and guest.")

    # Setup residue information
    mi1  =  Chem.AtomPDBResidueInfo()
    mi1.SetResidueName(resnames['hostresname'])
    mi1.SetResidueNumber(1)
    mi1.SetOccupancy(1.0)
    mi1.SetTempFactor(0.0)
    
    mi2  =  Chem.AtomPDBResidueInfo()
    mi2.SetResidueName(resnames['guestresname'])
    mi2.SetResidueNumber(2)
    mi2.SetOccupancy(1.0)
    mi2.SetTempFactor(0.0)
    
    # Add residue information to molecule
    count = 0
    for a in mol.GetAtoms():
        if count < CB_ATOM_NUM:
            a.SetMonomerInfo(mi1)
            count+=1
        elif count >= CB_ATOM_NUM:
            a.SetMonomerInfo(mi2)
            count+=1		
        
    return mol

def ammend_pdb_spacing(pdbfile):
    """ 
    Modifies pdb file in place to behave with xTB/AMBER
    """ 
    with open(pdbfile,"r") as f:
        lines = f.readlines()
    with open(pdbfile,"w") as f:
        # Write header
        f.write("comment1\n\n")
        for line in lines:
            species = line.split()[-1]
            f.write(line[:13] + species + "   " + line[13:])
    return
