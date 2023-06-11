# =============================================================================
# MAIN FOR SUPRAMOLECULAR FREE ENERGY OF BINDING CALCULATOR
# ============================================================================
# Author(s): Andy Thomson
# Configuration is changed via the flags and options in config.ini
# run() handles data in and out
# class methods handle module imports and running

import configparser
import importlib
import pandas as pd
from collections import OrderedDict
import subprocess as sp

from prog_ops import delete_and_create_dir

# Configuration manager class
class ConfigManager:
    """Class for managing the configuration of the program
    """

    def __init__(self, config_file_path):
        self.config_file_path = config_file_path
        self.config = configparser.ConfigParser(allow_no_value=True)
        self.config.read(config_file_path)

    def get_batchsize(self):
        """ Return batch size"""
        return int(self.config.get("INPUT", "batch"))

    def get_input_config(self):
        """ Return input condiguration that is global to all steps"""
        section = "INPUT"
        stepflags = self.config.get(section, "steps").split(",")
        try:
            df = pd.read_csv(self.config.get(section, "dataframefile"),sep="\t",names=['PubChemID','SMILES'])
        except:
            df = pd.read_pickle(self.config.get(section, "dataframefile"))
            print("Starting data:\n",df)
        df.index.name = 'MolId'
        return df, stepflags
    
    def get_molgen_config(self):
        """ Return molgen configuration"""
        section = "MOLGEN"
        method = self.config.get(section, "method")
        conformeroutdir = self.config.get(section, "conformeroutdir")
        return method, conformeroutdir
    
    def get_optimisation_config(self):
        """ Return optimisation configuration"""
        section = "OPTIMISATION"
        method = self.config.get(section, "method")
        inp = self.config.get(section, "inputfile")
        hostdir = self.config.get("INPUT", "hostdir")
        return method, inp, hostdir
    
    def get_docking_config(self):
        """ Return docking configuration"""
        section = "DOCKING"
        method = self.config.get(section, "method")
        inp = self.config.get(section, "inputfile")
        dockingoutdir = self.config.get(section, "dockingoutdir")
        return method, inp, dockingoutdir
    
    def get_final_outdir(self):
        """ Return final outdir"""
        return self.config.get("OUTPUT", "finaloutdir")

    def get_checks_flags(self):
        """ Return checks flags"""
        section = "CHECKS"
        conformercheck = self.config.getboolean(section, "conformercheck")
        bindingposecheck = self.config.getboolean(section, "bindingposecheck")
        exclusioncomplexcheck = self.config.getboolean(section, "exclusioncomplexcheck")
        return conformercheck, bindingposecheck, exclusioncomplexcheck

    
class ChemistrySimulator:
    """Class for managing the chemistry simulation
    """

    def __init__(self, config_file_path):
        self.config = ConfigManager(config_file_path)
        # Flag for whether the host has been optimised
        self.host_optimised = True

    # Run through each step with a batch of data at a time and update dataframe
    def run_batchwise(self, batchsize):
        """ Run through each step with a batch of data at a time and update dataframe"""

        # Retrieve data
        df, stepflags = self.config.get_input_config()

        # Check flags for checks
        conformercheck, bindingposecheck, exclusioncomplexcheck = self.config.get_checks_flags()

        # Create list of flags for each step
        steps = OrderedDict({"molgen":False,"optimisation":False,"docking":False,"complex_opt":False,"binding":False})
        for step in steps:
            if step in stepflags:
                steps[step] = True
                print("Running:",step)

        # MOL GENERATION
        # ================================================================
        molgen_choice, conformeroutdir = self.config.get_molgen_config()
        if steps["molgen"]:

            # Df starts with index: MolID; columns: PubChemID, SMILES
            batch = df["SMILES"].iloc[:batchsize]

            delete_and_create_dir(f"{conformeroutdir}/confs")

            # Generate conformers
            conformers = []
            for i in range(batchsize):
                conformer = self.model_gen(molgen_choice,
                                           batch[i],
                                           df.index[i],
                                           outdir=f"{conformeroutdir}/confs"
                                           )
                conformers.append({"MolID":df.index[i],"guestmol":conformer})
        
            conformers = pd.DataFrame(conformers)
            conformers.set_index("MolID",inplace=True)

            df = pd.concat([df,conformers],axis=1)
            df.to_pickle(conformeroutdir + "/conformerdf.pkl")

        # OPTIMISATION
        # ================================================================
        opt_choice, inp, hostdir = self.config.get_optimisation_config()
        if steps["optimisation"]:

            # Df starts with index: MolID; columns: PubChemID, SMILES, guestmol
            batch = df["guestmol"].iloc[:batchsize]

            delete_and_create_dir(f"{conformeroutdir}/opt")

            # Optimise guest
            optguests = []
            for i in range(batchsize):
                if batch[i] != "InvalidSMILES":
                    optguest = self.optimise(opt_choice,
                                             batch[i],
                                             df.index[i],
                                             inp,
                                             outdir=f"{conformeroutdir}/opt"
                                             )
                    optguests.append({"MolId":df.index[i],"guestmol":optguest})
                else:
                    optguests.append({"MolId":df.index[i],"guestmol":"InvalidSMILES"})

            optguests = pd.DataFrame(optguests)
            optguests.set_index("MolId",inplace=True)
            
            df.update(optguests)
            df.to_pickle(conformeroutdir + "/conformerdf.pkl")

        # DOCKING
        # ================================================================
        docking_choice, inp, dockingoutdir = self.config.get_docking_config()
        if steps["docking"]:

            hostfile = f"{hostdir}/xtbopt.pdb"

            # Df starts with index: MolID; columns: PubChemID, SMILES, guestmol
            batch = df["guestmol"].iloc[:batchsize]

            delete_and_create_dir(f"{dockingoutdir}/docked")

            # Dock guest to host
            complexes = []
            for i in range(batchsize):
                if batch[i] != "InvalidSMILES":
                    compl = self.docking(docking_choice,
                                         batch[i],
                                         df.index[i],
                                         hostfile,
                                         None,
                                         outdir=f"{dockingoutdir}/docked"
                                         )
                    complexes.append({"MolId":df.index[i],"dockedmol":compl})
                else:
                    complexes.append({"MolId":df.index[i],"dockedmol":"InvalidSMILES"})

            complexes = pd.DataFrame(complexes)
            complexes.set_index("MolId",inplace=True)

            df = pd.concat([df,complexes],axis=1)
            df.to_pickle(dockingoutdir + "/dockingdf.pkl")

        # COMPLEX OPTIMISATION
        # ================================================================
        opt_choice, inp, hostdir = self.config.get_optimisation_config()
        if steps["docking"]:

            # Df starts with index: MolID; columns: PubChemID, SMILES, guestmol, dockedmol
            batch = df["dockedmol"].iloc[:batchsize]

            delete_and_create_dir(f"{dockingoutdir}/complexopt")

            # Optimise complex
            optcomplexes = []
            for i in range(batchsize):
                if batch[i] != "InvalidSMILES":
                    optcomplex = self.optimise(opt_choice,
                                               batch[i],
                                               df.index[i],
                                               inp,
                                               outdir=f"{dockingoutdir}/complexopt"
                                                )
                    optcomplexes.append({"MolId":df.index[i],"dockedmol":optcomplex})
                else:
                    optcomplexes.append({"MolId":df.index[i],"dockedmol":"InvalidSMILES"})

            optcomplexes = pd.DataFrame(optcomplexes)
            optcomplexes.set_index("MolId",inplace=True)

            df.update(optcomplexes)
            df.to_pickle(dockingoutdir + "/dockingdf.pkl")

        # EXCLUSION COMPLEX
        # ================================================================
        if exclusioncomplexcheck:

            batch = df["dockedmol"].iloc[:batchsize]

            # Check for exclusion complex
            exclusioncomplexes = []
            for i in range(batchsize):
                if batch[i] != "InvalidSMILES":
                    exo = {}
                    isExo, centroid_diff, cavity_atoms = self.check_exclusion_complex(batch[i])
                    exo["isExo"] = isExo
                    exo["centroid_diff"] = centroid_diff
                    exo["cavity_atoms"] = cavity_atoms
                    exclusioncomplexes.append({"MolId":df.index[i]}.update(exo))
                else:
                    exclusioncomplexes.append({"MolId":df.index[i],"isExo":"InvalidSMILES"})

            exclusioncomplexes = pd.DataFrame(exclusioncomplexes)
            exclusioncomplexes.set_index("MolId",inplace=True)

            df = pd.concat([df,exclusioncomplexes],axis=1)
            df.to_pickle(dockingoutdir + "/dockingdf.pkl")

        # BINDING ENERGY
        # ================================================================
        binding_choice = opt_choice
        finaloutdir = self.config.get_final_outdir()
        if steps["binding"]:

            hostoutfile = f"{hostdir}/opt.out"

            # Df starts with index: MolID; columns: PubChemID, SMILES, guestmol, dockedmol
            batch = df["dockedmol"].iloc[:batchsize]
            batchguest = df["guestmol"].iloc[:batchsize]

            outfilealreadyexist = True
            # Check dockingoutdir and conformeroutdir outfiles
            if open(f"{dockingoutdir}/dockingdf.pkl","rb") == None or open(f"{conformeroutdir}/conformerdf.pkl","rb") == None:
                complexmolfinals = []
                guestmolfinals = []
                outfilealreadyexist = False
            
            # Calculate binding energy
            bindingenergies = []
            for i in range(batchsize):
                if batch[i] != "InvalidSMILES":
                    if outfilealreadyexist:
                        complexoutfile = f"{dockingoutdir}/{df.index[i]}_opt.out"
                        guestoutfile = f"{conformeroutdir}/{df.index[i]}_opt.out"
                        bindingenergy = self.binding(binding_choice,
                                                    batch[i],
                                                    batchguest[i],
                                                    df.index[i],
                                                    hostoutfile,
                                                    inp,
                                                    complexoutfile=complexoutfile,
                                                    guestoutfile=guestoutfile,
                                                    outdir=finaloutdir
                                                    )
                        bindingenergies.append({"MolId":df.index[i]}.update(bindingenergy))
                    else:
                        complexmolfinal, guestmolfinal, bindingenergy = self.binding(binding_choice,
                                                    batch[i],
                                                    batchguest[i],
                                                    df.index[i],
                                                    hostoutfile,
                                                    inp,
                                                    outdir=finaloutdir
                                                    )
                        bindingenergies.append({"MolId":df.index[i]}.update(bindingenergy))
                        complexmolfinals.append({"MolId":df.index[i],"dockedmol":complexmolfinal})
                        guestmolfinals.append({"MolId":df.index[i],"guestmol":guestmolfinal})
                else:
                    bindingenergies.append({"MolId":df.index[i],"bindingenergy":"InvalidSMILES"})

            if not outfilealreadyexist:
                complexmolfinals = pd.DataFrame(complexmolfinals)
                complexmolfinals.set_index("MolId",inplace=True)
                guestmolfinals = pd.DataFrame(guestmolfinals)
                guestmolfinals.set_index("MolId",inplace=True)
                df.update(complexmolfinals)
                df.update(guestmolfinals)

            bindingenergies = pd.DataFrame(bindingenergies)
            bindingenergies.set_index("MolId",inplace=True)

            df = pd.concat([df,bindingenergies],axis=1)
            df.to_pickle(finaloutdir + "/finaldf.pkl")

        # DATAFILTER
        # ================================================================

        return
            
    # Run through a single piece of data and return a single piece of data (for integration into MCTS)
    def run_piecewise(self):
        #file = self.config.get_input_config() 

        # Mol generation

        # Optimisation

        # Docking

        # Binding energy

        pass

    def run(self):
        batch = self.config.get_batchsize()
        if batch != 1:
            self.run_batchwise(batch)
        else:
            self.run_piecewise()

    # SMILES to 3D structure
    def model_gen(self,method,smiles,molId,outdir):
        """ Steps: convert smiles to mol object
        Embed mol object to get 3D structure
        """
        module = importlib.import_module(f"molgen_{method}")
        mol = module.mol_gen(smiles,molId)

        # Files created:
        #f"{molId}_rdkit.pdb"
        #f"{molId}_reduced.pdb"

        # Move files
        sp.run["mv",f"{molId}_{method}.pdb",f"{outdir}/{molId}_{method}.pdb"]
        sp.run["mv",f"{molId}_reduced.pdb",f"{outdir}/{molId}_reduced.pdb"]

        return mol

    # Guest and host optimisation
    def optimise(self,method,mol,molId,inp,outdir):
        """Steps: optimise guest with chosen ESM (electronic structure method) options
        Optimise host with chosen ESM options if host_optimsed is False
        """
        module = importlib.import_module(f"opt_{method}")
        guestmol = module.opt(mol,molId,inp)

        # Files created:
        #f"{molId}_opt.out"
        #f"{molId}_opt.pdb"

        # Move files
        sp.run["mv",f"{molId}_opt.out",f"{outdir}/{molId}_opt.out"]
        sp.run["mv",f"{molId}_opt.pdb",f"{outdir}/{molId}_opt.pdb"]
            
        return guestmol

    def docking(self,method,mol,molId,hostfile,inp,outdir):
        """Steps: dock guest to host
        # optimise complex
        # optmiise guest if config is different
        """
        module = importlib.import_module(f"dock_{method}")
        complexmol = module.dock(mol,molId,hostfile,inp)

        # Files created:
        #f"{molId}_dock.out"
        #f"{molId}_complex.pdb"

        # Move files
        sp.run["mv",f"{molId}_dock.out",f"{outdir}/{molId}_dock.out"]
        sp.run["mv",f"{molId}_complex.pdb",f"{outdir}/{molId}_complex.pdb"]

        return complexmol

    def binding_energy(self,method,complexmol,guestmol,molId,hostoutfile,inp,**kwargs):
        """Steps: calculate binding energy
        Get energy breakdown (kwargs could indicate which breakdown to get)
        """
        module = importlib.import_module(f"bind_{method}")

        # Returns dictionary with breakdown of binding energy
        if kwargs["complexoutfile"] and kwargs["guestoutfile"]:
            bindingdict = module.binding_energy(complexmol,guestmol,molId,hostoutfile,inp,complexoutfile=kwargs["complexoutfile"],guestoutfile=kwargs["guestoutfile"])
        else:
            bindingdict = module.binding_energy(complexmol,guestmol,molId,hostoutfile,inp,outdir=kwargs["outdir"])

        return bindingdict
    
    def check_exclusion_complex(self,complexmol):
        """Steps: check if complex is excluded
        """
        module = importlib.import_module("exclusion.py")

        isExo, centroid_diff, cavity_atoms = module.exclusion_check(complexmol)

        return isExo, centroid_diff, cavity_atoms

    def thermo(self):
        """Steps: extensible method for calculating entropic contribution to free energy of binding
        DATA CHECKS:
        Check for convergence
        Return thermodynamic properties
        """
        pass
        
if __name__ == "__main__":
    config_file_path = "config.ini"
    simulator = ChemistrySimulator(config_file_path)
    simulator.run()
    # print(simulator.config.get_input_config())
    # print(simulator.config.get("MOLGEN", "method"))
    # print(simulator.config.get("OPTIMISATION", "method"))
    # print(simulator.config.get("DOCKING", "method"))
    # print(simulator.config.get("BINDING_ENERGY", "method"))
    # print(simulator.config.get("THERMO", "method"))
    # print(simulator.config.get("INPUT", "file"))
    # print(simulator.config.get("OUTPUT", "file"))
