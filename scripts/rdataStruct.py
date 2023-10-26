################################################################################################
# @info Classes for the ROOT input/output files from the various devices                       #
# @creation date 23/04/19                                                                      #
# @edit     date 23/08/06                                                                      #
#                                               more text here eventually                      #
################################################################################################
from logger import create_logger
from datetime import datetime
import numpy as np
import ROOT
import tqdm





######################################################
######################################################
######################################################
class rdataStruct_OPT():
    # rdataStruct_FERS logger
    logging = create_logger("rdataStruct_OPT")

    def __init__(self, fname, mode="RECREATE", setVars = None):
        ## Internal
        # mode state
        self.mode = mode
        # filename
        self.ROOTfilename = fname
    
        ######################################################
        ################
        ### TTree: OPT
        ## Variable definitions
        #self.FERS = None
        self.OPT_ibunch = np.array([0], dtype=np.uint32)
        self.OPT_ievt = np.array([0], dtype=np.uint32)
        self.OPT_idet = np.array([0], dtype=np.uint32)
        self.OPT_idir = np.array([0], dtype=np.uint32)
        self.OPT_voptPar    = np.zeros(3, dtype=np.double)
        self.OPT_voptParErr = np.zeros(3, dtype=np.double)

        self.OPT_nametypes = [
            ("bunch",           self.OPT_ibunch,        "bunch/i",              "Bunch number"),
            ("evt",             self.OPT_ievt,          "evt/i",                "Event number"),
            ("det",             self.OPT_idet,          "det/i",                "Detector ID [0-1]"),
            ("dir",             self.OPT_idir,          "dir/i",                "Direction [0 for X, 1 for Y]"),
            #
            ("optPar",          self.OPT_voptPar,       "optPar[3]/D",          "Profile fit parameters"),
            ("optParErr",       self.OPT_voptParErr,    "optParErr[3]/D",       "Profile fit parameter errors")
        ]
        self.OPT_fill_warnings = {item[0]:False for item in self.OPT_nametypes}  
        ######################################################
    

        # Open ROOT output file
        self.ROOTfile = ROOT.TFile.Open(fname, mode)
        (self.logging).debug(f"Opening {fname} in {mode} mode")
        if mode == "RECREATE":
            (self.logging).debug("Attaching branches to variables.")
            # Generate the OPT TTree
            self.setupTTree_OPT()
        else:
            # OPT
            self.OPT = self.ROOTfile.OPT
            # Set branch addresses
            self.OPT_bindBranches(self.OPT)
    
    ## Functions
    # Bind branches of the OPT TTree to variables
    def OPT_bindBranches(self, tree):
        for item in self.OPT_nametypes:
            tree.SetBranchAddress(item[0], item[1])
    # Fill tree with data
    def OPT_fill(self, **kwargs):
        if (self.OPT) is None:
            raise Exception("Tree 'data' is not initialized")
        for branch in self.OPT_nametypes:
            try:
                # Filling is based on length. If the length of the array is 1 then fill the content of the array
                # otherwise fill the array itself
                if len(branch[1]) == 1:
                    branch[1][0] = kwargs[branch[0]]
                else:
                    np.copyto(branch[1], kwargs[branch[0]])
            except KeyError:
                if not self.OPT_fill_warnings[branch[0]]:
                    (self.logging).warning(f"Leaf {branch[0]} not found in kwargs (further warnings suppressed)")
                    self.OPT_fill_warnings[branch[0]] = True
        # Fill the tree
        (self.OPT).Fill()
    # Get entry i-th of the OPT TTree (for read mode)
    def OPT_getEntry(self, i):
        """
        Description 
        
        Parameters
        ----------
            i (int): Entry number 

        Returns:
        
            OPT_ibunch : 0
            OPT_ievt : 1
            OPT_idet : 2
            OPT_idir : 3
            OPT_voptPar : 4
            OPT_voptParErr : 5
        """
        (self.OPT).GetEntry(i)
        return [(self.OPT_ibunch)[0], (self.OPT_ievt)[0], (self.OPT_idet)[0], (self.OPT_idir)[0], self.OPT_voptPar, self.OPT_voptParErr]
    # Setup the OPT TTree
    def setupTTree_OPT(self):
        # Define the TTree
        self.OPT = ROOT.TTree("OPT", "Digitization opt analysis data")
        # Create branches
        for item in (self.OPT_nametypes):
            (self.OPT).Branch(item[0], item[1], item[2])
        # Set leaves description
        self.TREE_SetLeavesDescriptions(self.OPT, self.OPT_nametypes)

    
    
    ######################################################
    ######################################################
    ## Utilities
    # Set description (better readibility)
    def TREE_SetLeavesDescriptions(self, tree, tree_nametypes):
        for entry in tree_nametypes:
            tree.GetBranch(entry[0]).SetTitle(entry[3])

    # Write all the content of the output root file and clear datastructures for the next one
    def closeROOT(self):
        if self.mode != "READ":
            (self.OPT).Write()
        # Single file
        if type(self.ROOTfilename) is not list:
            (self.ROOTfile).Close()






######################################################
######################################################
######################################################

# Test methods
if __name__ == "__main__":
    # rdataStruct_FERS logger
    logging = create_logger("test_rdataStruct", 10)
    
    
    # Test the MERGE class
    logging.info("Testing rdataStructMERGE class by generating a testMERGE.root in this dir...")
    roFile = rdataStruct_OPT("testBuddyPLAIN.root")
    logging.info("Fill 1 without all kwargs")
    roFile.OPT_fill(bunch = 0, evt = 0)
    logging.info("Fill 2 without all kwargs")
    roFile.OPT_fill(bunch = 0, evt = 1)
    logging.info("Fill 2 witg all kwargs")
    roFile.OPT_fill(bunch = 0, evt = 2, det=0, dir=0, optPar=np.array([0,1,2]), optParErr = np.array([0,1,2]))
    roFile.closeROOT()
    logging.info("File closed.")