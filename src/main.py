# =============================================================================
# MAIN FOR SUPRAMOLECULAR FREE ENERGY OF BINDING CALCULATOR
# ============================================================================
# Author(s): Andy Thomson
# Configuration is changed via the flags and options in config.ini

import configparser
import importlib
import pandas as pd
from collections import OrderedDict
import subprocess as sp

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
    
class ChemistrySimulator:
    """Class for managing the chemistry simulation
    """

    def __init__(self, config_file_path):
        self.config = ConfigManager(config_file_path)
        # Flag for whether the host has been optimised
        self.host_optimised = False

    # Run through each step with a batch of data at a time and update dataframe
    def run_batchwise(self, batchsize):
        """ Run through each step with a batch of data at a time and update dataframe"""

        # Retrieve data
        df, stepflags = self.config.get_input_config()

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

            # Generate conformers
            conformers = []
            for i in range(batchsize):
                conformers.append({"MolID":df.index[i],"guestmol":self.model_gen(molgen_choice,batch[i],df.index[i],outdir=conformeroutdir)})
        
            conformers = pd.DataFrame(conformers)
            conformers.set_index("MolID",inplace=True)

            df = pd.concat([df,conformers],axis=1)
            df.to_pickle(conformeroutdir + "/conformerdf.pkl")

        # OPTIMISATION
        # ================================================================
        opt_choice, inp, hostdir = self.config.get_optimisation_config()
        if steps["optimisation"]:        

            # Check if host has been optimised
            # This is xTB specific, will have to change
            if open(f"{hostdir}/opt.out", "r").read() == "True":
                self.host_optimised = True

            # If host has not been optimised, optimise it
            if not self.host_optimised:
                pass

            # Df starts with index: MolID; columns: PubChemID, SMILES, guestmol
            batch = df["guestmol"].iloc[:batchsize]

            # Optimise guest
            optguests = []
            for i in range(batchsize):
                if batch[i] != "InvalidSMILES":
                    optguests.append({"MolId":df.index[i],"guestmol":self.optimise(opt_choice,batch[i],df.index[i],inp,outdir=conformeroutdir)})
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

            # Dock guest to host
            complexes = []
            for i in range(batchsize):
                if batch[i] != "InvalidSMILES":
                    complexes.append({"MolId":df.index[i],"dockedmol":self.dock(docking_choice,batch[i],df.index[i],hostfile,None,outdir=dockingoutdir)})
                else:
                    complexes.append({"MolId":df.index[i],"dockedmol":"InvalidSMILES"})

            complexes = pd.DataFrame(complexes)
            complexes.set_index("MolId",inplace=True)

            df = pd.concat([df,complexes],axis=1)
            df.to_pickle(dockingoutdir + "/dockingdf.pkl")

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
    def model_gen(self,method,smiles,molId,**kwargs):
        """ Steps: convert smiles to mol object
        Embed mol object to get 3D structure
        DATA CHECKS:
        Check for different conformers
        Ensure that the structure is hydrogenated
        KWARGS:
        Choose checks and outdir for molecule generated
        Return mol object
        """
        module = importlib.import_module(f"molgen_{method}")
        if not kwargs["outdir"]:
            kwargs["outdir"] = None
        mol = module.single_smiles(smiles,molId,outdir=kwargs["outdir"])

        return mol

    # Guest and host optimisation
    def optimise(self,method,mol,molId,inp,**kwargs):
        """Steps: optimise guest with chosen ESM (electronic structure method) options
        Optimise host with chosen ESM options if host_optimsed is False
        DATA CHECKS:
        Check for convergence
        KWARGS:
        Choose checks, outdir if they exist
        Return mol object
        """
        module = importlib.import_module(f"opt_{method}")

        if not kwargs["outdir"]:
            outfile = f"{molId}_opt.out"
            optfile = f"{molId}_opt.pdb"

        
        else:
            outfile = kwargs["outdir"] + f"/{molId}_opt.out"
            optfile = kwargs["outdir"] + f"/{molId}_opt.pdb"
            
        guestmol = module.opt(mol,molId,outfile,optfile,inp)
            
        # Clean up files
        if not kwargs["outdir"]:
            sp.run(["rm",outfile,optfile])

        return guestmol

    def docking(self,method,mol,molId,hostfile,inp,**kwargs):
        """Steps: dock guest to host
        # optimise complex
        # optmiise guest if config is different
        DATA CHECKS:
        Check for exclusion complexes
        Check for multiple docking poses
        Check convergence
        Return mol object
        """
        module = importlib.import_module(f"dock_{method}")

        if not kwargs["outdir"]:
            dockoutfile = f"{molId}_dock.out"
            posefile = f"{molId}_complex.pdb"

        else:
            dockoutfile = kwargs["outdir"] + f"/{molId}_dock.out"
            posefile = kwargs["outdir"] + f"/{molId}_complex.pdb"

        # This is xTB specific, will have to change
        complexmol = module.dock(mol,molId,hostfile,dockoutfile,posefile,inp)

        # Clean up
        if not kwargs["outdir"]:
            sp.run(["rm",dockoutfile,posefile])

        return complexmol

    def binding_energy(self):
        """Steps: calculate binding energy
        Get energy breakdown
        DATA CHECKS:
        Flag if binding energy is very high or very low
        Return binding energy"""
        pass

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
