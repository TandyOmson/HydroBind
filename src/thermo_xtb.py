""" Module for calculating normal mode frequencies using xTB
xTB uses SSHRO to calculate normal modes, entropy is extracted
There should be no imaginary frequencies (negative Hessian eigenvalues)
"""

from rdkit import Chem
import subprocess as sp
from mol_ops import ammend_pdb_spacing

def read_thermo(outfile):
    """ reads in the thermochemistry output from xTB """

    en_dict = {}
    # Get final energy breakdown from end of file
    gen = (i.split() for i in reversed(open(outfile,"r").readlines()))

    for i in gen:
        if " ".join(i[:3]) == ":: G(RRHO) contrib.":
            en_dict["G_rrho"]  = float(i[3]) * 627.509
        if " ".join(i[:1]) == "TR":
            en_dict["H_tr"] = float(i[2]) * 10E-3
            en_dict["TS_tr"] = float(i[4]) * 10E-3 * 298.15
        if " ".join(i[:1]) == "INT":
            en_dict["H_int"] = float(i[2]) * 10E-3
            en_dict["TS_int"] = float(i[4]) * 10E-3 * 298.15
        if " ".join(i[:1]) == "ROT":
            en_dict["H_rot"] = float(i[2]) * 10E-3
            en_dict["TS_rot"] = float(i[4]) * 10E-3 * 298.15
        if " ".join(i[:2]) == "298.15 VIB":
            en_dict["H_vib"] = float(i[3]) * 10E-3
            en_dict["TS_vib"] = float(i[5]) * 10E-3 * 298.15
            break

    return en_dict

def thermo(complexmol,guestmol,molId,hostfile):
    """ optimises a given molecule using xTB """
    complexinfile = f"{molId}.pdb"
    guestinfile = f"{molId}_guest.pdb"
    complexoutfile = f"{molId}_hess.out"
    guestoutfile = f"{molId}_guest_hess.out"

    Chem.MolToPDBFile(complexmol,f"{complexinfile}")
    # MolToPDBFile formats differently for complexes and single molecules, this is a temporary fix
    #if len(complexmol.GetAtoms()) > 125:
    #    ammend_pdb_spacing(f"{complexinfile}")

    Chem.MolToPDBFile(guestmol,f"{guestinfile}")

    sp.run(["xtb",f"{guestinfile}","--hess","--alpb","water"],stdout=open(guestoutfile,"w"))
    sp.run(["xtb",f"{complexinfile}","--hess","--alpb","water"],stdout=open(complexoutfile,"w"))

    # Read in thermochemistry
    complex_thermo = read_thermo(complexoutfile)
    guest_thermo = read_thermo(guestoutfile)
    host_thermo = read_thermo(hostfile)

    thermo_dict = {
                    key: complex_thermo[key] - guest_thermo[key] - host_thermo[key]
                    for key in complex_thermo.keys()
                    }
    
    thermo_dict["MolId"] = str(molId)

    # Clean up
    sp.run(["rm",f"{complexinfile}"])
    sp.run(["rm",f"{guestinfile}"])
    sp.run(["rm",f"{complexoutfile}"])
    sp.run(["rm",f"{guestoutfile}"])
    sp.run(["rm","hessian"])
    sp.run(["rm","g98.out"])
    sp.run(["rm","vibspectrum"])
    sp.run(["rm","wbo"])
    sp.run(["rm","xtbtopo.mol"])
    sp.run(["rm","xtbrestart"])
    sp.run(["rm","charges"])

    return thermo_dict