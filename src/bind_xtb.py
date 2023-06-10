# Chemistry simulator module. 
# Uses xTB to generate free energy of binding for a given complex

from rdkit import Chem
import subprocess as sp
import pandas as pd

from mol_ops import ammend_pdb_spacing

def get_binding_from_outfile(outfile,molId):
    """ Gets binding energy from xTB output file
    Returns a dictionary of binding energies
    Could add kwargs to specify which energies to return using a dictionary of bools
    """
    
    en_dict = {}
    # Get final energy breakdown from end of file
    gen = (i.split() for i in reversed(open(outfile,"r").readlines()))

    # Get molecule name
    en_dict['mol'] = f"D{str(molId)}"
    
    for i in gen:
        if " ".join(i[:3]) == ":: repulsion energy":
            en_dict["repulsion"]  = float(i[3])

        if " ".join(i[:3]) == ":: -> Gshift":
            en_dict["Gshift"] = float(i[3])
            
        if " ".join(i[:3]) == ":: -> Ghb":
            en_dict["Ghb"] = float(i[3])
            
        if " ".join(i[:3]) == ":: -> Gsasa":
            en_dict["Gsasa"] = float(i[3])
                
        if " ".join(i[:3]) == ":: -> Gelec":
            en_dict["Gelec"] = float(i[3])

        if " ".join(i[:3]) == ":: -> Gsolv":
            en_dict["Gsolv"] = float(i[3])

        if " ".join(i[:3]) == ":: -> dispersion":
            en_dict["disp"] = float(i[3])
        
        if " ".join(i[:4]) == ":: -> anisotropic XC":
            en_dict["aniso_XC"] = float(i[4])

        if " ".join(i[:4]) == ":: -> anisotropic ES":
            en_dict["aniso_ES"] = float(i[4])
            
        if " ".join(i[:4]) == ":: -> isotropic ES":
            en_dict["iso_ES"] = float(i[4])
            
        if " ".join(i[:3]) == ":: SCC energy":
            en_dict["SCC"]= float(i[3])

        if " ".join(i[:3]) == ":: total energy":
            en_dict["toten"] = float(i[3])
            break
                        
    for i,j in en_dict.items():
        if not isinstance(j,str):
            en_dict[i] = en_dict[i]*627.509

    return en_dict

def bind(complexmol,molId,hostoutfile,xtbinp,**kwargs):
    """ calculates binding affinity of a given complex using xTB
    kwargs include both the optional outdir and optional guestoutfile (output from guest optimisation)
    if notguestoutfile, then the guestpdb is given.
    Returns a dictionary for each molecule of all relevant binding information
    """
    # Gets a breakdown of xTB energies, returns it as a dictionary. Requires outfiles from complex, guest and host SCFs
    
    # xTB on complex
    if not kwargs["guestoutfile"] or kwargs["complexoutfile"]:
        complexoutfile = f"{molId}_comp_opt.out"
        Chem.MolToPDBFile(complexmol,f"{molId}_comp.pdb")
        ammend_pdb_spacing(f"{molId}_comp.pdb")
        sp.run(["xtb","--input",f"{xtbinp}",f"{molId}_comp.pdb","--opt","--alpb","water"],stdout=open(complexoutfile,"w"))
        # Get final mol objects for complex
        complexmol_final = Chem.MolFromPDBFile("xtbopt.pdb",removeHs=False,sanitize=False)

        guestoutfile = f"{molId}_opt.out"
        sp.run(["xtb","--input",f"{xtbinp}",f"{kwargs['guestpdb']}","--opt","--alpb","water"],stdout=open(guestoutfile,"w"))
        # Get final mol object for guest
        guestmol_final = Chem.MolFromPDBFile("xtbopt.pdb",removeHs=False,sanitize=False)

    else:
        complexoutfile = kwargs["complexoutfile"]
        guestoutfile = kwargs["guestoutfile"]

    # Calculate binding affinities
    guest_en_dict = get_binding_from_outfile(guestoutfile,molId)
    host_en_dict = get_binding_from_outfile(hostoutfile,molId)
    complex_en_dict = get_binding_from_outfile(complexoutfile,molId)

    binding_dict = {
                    key: complex_en_dict[key] - guest_en_dict[key] - host_en_dict[key]
                    if not isinstance(complex_en_dict[key], str)
                    and not isinstance(guest_en_dict[key], str)
                    and not isinstance(host_en_dict[key],str)
                    else complex_en_dict[key] for key in complex_en_dict
                    }

    # Move outfiles to outdir if given
    if not kwargs["guestoutfile"] or kwargs["complexoutfile"]:

        if kwargs.get('outdir'):
            sp.run(["mv",complexoutfile,kwargs['outdir']])
            sp.run(["cp",guestoutfile,kwargs['outdir']])

        # Clean up
        sp.run(["rm","xtbopt.pdb"])
        sp.run(["rm",f"{molId}_comp.pdb"])
        sp.run(["rm","xtbtopo.mol"])
        sp.run(["rm","wbo"])
        sp.run(["rm","xtbrestart"])
        sp.run(["rm","xtbopt.log"])
        sp.run(["rm",".xtboptok"])
        sp.run(["rm","charges"])

    return complexmol_final, guestmol_final, binding_dict

# To move to test
if __name__ == "__main__":
    # Grab current data from DataFrame
    df = pd.read_pickle("/home/spine/DProjects/DhydroEn/results/Dhydrocarbons/all_mols_3.pkl")

    # xTB input file
    xtbinp = "/home/spine/DProjects/DhydroEn/data/xtbinp.inp"
    # Host optimisation output file
    hostoutfile = "/home/spine/DProjects/DhydroEn/data/DCB7/opt.out"

    # Get complex mol and molId
    complexmol = df["docked_mol"].values
    molId = df.index.values

    complexmoldict = {}
    guestmoldict = {}

    binding_df = pd.DataFrame()
    k = 0
    while k != 1:
        for i,j in zip(complexmol,molId):
            complexmol_final, guestmol_final, binding_dict = bind(
                          i,
                          j,
                          hostoutfile,
                          xtbinp,
                          outdir="/home/spine/DProjects/DhydroEn/results/DfinalHydrocarbons",
                          guestoutfile=f"/home/spine/DProjects/DhydroEn/results/DfinalHydrocarbons/{j}_opt.out"
                          )

            complexmoldict[j] = complexmol_final
            guestmoldict[j] = guestmol_final
            bindingdict = pd.DataFrame.from_dict(binding_dict,orient="index")
            print(bindingdict)
            binding_df = pd.concat([binding_df,bindingdict],axis=1)

        k += 1

    complexdf = pd.DataFrame.from_dict(complexmoldict,orient="index",columns=["docked_mol"])
    guestdf = pd.DataFrame.from_dict(guestmoldict,orient="index",columns=["mol"])
    # Update docked_mol and mol columns with optimised mol objects
    df.update(complexdf)
    df.update(guestdf)

    # Update binding energies
    df = pd.concat([df,binding_df.T],axis=1)
    print(df)

    df.to_pickle("/home/spine/DProjects/DhydroEn/results/Dhydrocarbons/all_mols_4.pkl")
