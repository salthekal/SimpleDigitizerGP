#################################################################################################
# @info description                                                                 #
# @date 23/10/26                                                                                #
#                                                                                               #

#################################################################################################
from multipledispatch import dispatch
from logger import create_logger
from tqdm import tqdm
import numpy as np
import ROOT
import os


class readFromMc():
    # readFromMc class logger
    logging = create_logger("readFromMc")
    
    # Electron charge in coulomb
    eCharge = 1.602176634e-19
    
    
    def __init__(self, rFname: str, bunchParNb: int, cce: float, avgPairEn: float):
        """
        DD
        
        Parameters
        ----------
            rFname (str) : Filename of the MC ROOT file containing the data
            bunchParNb (int) : number of particles in a bunch
            cce (float) : Charge collection efficiency of the sensor for the bare geometrical projection of the dep chg. to proj chg.
            avgPairEn (float) : Average energy required to create an electron-hole pair in the sensor (in eV)
        
        Returns
        -------
            None
        """
        
        # Initialize the internal class variables with the given ext ones.
        self.rFname = rFname
        self.bunchParNb = bunchParNb
        self.cce = cce
        self.avgPairEn = avgPairEn
        
        # Checks validity of the parameters
        if not (cce>=0 and cce<=1): raise Exception("CCE must be between 0 and 1")
        
        # Enable ROOT implicit multithreading
        #ROOT.EnableImplicitMT(16)
    


    # Takes the map of energy deposited in a sensor and returns the map of charge collected at the strips (hor./vert.)
    def getChgDepProjProfiles(self, edepHist: ROOT.TH2D) -> dict:
        """
        Takes the map of energy deposited in a sensor and returns the map of charge collected at the strips (hor./vert.)
        
        Parameters
        ----------
            edepHist (ROOT.TH2D) : histogram with the map of energy deposition in a sensor
            cce (float) : Charge collection efficiency of the sensor for the bare geometrical projection of the dep chg. to proj chg.
            avgPairCreationEnergy_eV (float) : Average energy required to create an electron-hole pair in the sensor (in eV)
            
        Returns
        -------
            'chgDepProfile': (chgDepProfileX.Clone(), chgDepProfileY.Clone()),
            'chgProjProfile': (chgProjProfileX.Clone(), chgProjProfileY.Clone())
        """

        # Takes the map of energy deposition of det 'det' and returns the charge collected at the strip
        def projectCharge(chgDepProfile: ROOT.TH1D) -> ROOT.TH1D:
            """
            Takes the map of energy deposition of det 'det' and returns the charge collected at the strip
            
            Parameters
            ----------
                chgDepProfile (ROOT.TH1D) : histogram with the map of the charge deposited in the strips of the sensor 
            Returns
            -------
                chgProjProfile (ROOT.TH1D) : histogram with the map of the charge projected in the strips of the sensor 
            """
            
            # TODO: understand how the errors from the addition of the poisson in line 282 propagates with the errors in the chgProjProfile
            chgProjProfile = chgDepProfile.Clone()
            # Applies smearing of the deposited charge by adding some value draw from a Poissonian distribution also propagates the error on the deposited charge
            for i in range(chgProjProfile.GetNbinsX()):
                chgDep = chgProjProfile.GetBinContent(i)
                #chgDepError = chgProjProfile.GetBinError(i)
                chgDep += np.random.poisson(chgDep)
                #chgDepError = np.random.poisson(chgDepError)
                chgProjProfile.SetBinContent(i, chgDep)
                #chgDepProfile.SetBinError(i, chgDepError)   

            # Project to final charge at strips by multiplying for the CCE
            chgProjProfile.Scale(self.cce)
            # Get input title
            in_title = chgDepProfile.GetTitle()

            # Extract sensor name
            name_parts = in_title.split()
            sensor = name_parts[4]
        
            # Construct output title with same sensor name
            out_title = " ".join([name_parts[0], "projected", name_parts[2], name_parts[3], sensor, name_parts[5]])

            # Set output title
            chgProjProfile.SetTitle(out_title) 
            chgProjProfile.GetYaxis().SetTitle("charge created [C]")
            chgProjProfile.GetXaxis().CenterTitle()
            chgProjProfile.GetYaxis().CenterTitle()
            chgProjProfile.SetName(chgProjProfile.GetName().replace("Dep", "Proj"))
            return chgProjProfile


        # Get horizontal and vertical profiles of energy depositions in the sensor
        edepProfileX, edepProfileY = edepHist.ProjectionX(), edepHist.ProjectionY()
            
        # This maps the profile of edeps to a profile of charge deposited (in C)
        # The 1e3 accounts for the energy in the histo being in keV
        
        chgDepProfileX, chgDepProfileY = edepProfileX.Clone(), edepProfileY.Clone()
        chgDepProfileX.Scale(1e3 * self.eCharge / self.avgPairEn), chgDepProfileY.Scale(1e3 * self.eCharge / self.avgPairEn)

        # Change the titles in case one is interested to dump the histos here
        chgDepProfileX.SetTitle(chgDepProfileX.GetTitle().replace("Energy deposition", "Charge deposited"))
        chgDepProfileY.SetTitle(chgDepProfileY.GetTitle().replace("Energy deposition", "Charge deposited"))
        chgDepProfileX.GetYaxis().SetTitle("charge deposited [C]")
        chgDepProfileY.GetYaxis().SetTitle("charge deposited [C]")
        chgDepProfileX.GetXaxis().CenterTitle()
        chgDepProfileX.GetYaxis().CenterTitle()
        chgDepProfileY.GetXaxis().CenterTitle()
        chgDepProfileY.GetYaxis().CenterTitle()
        chgDepProfileX.SetName("chgDepProfileX")
        chgDepProfileY.SetName("chgDepProfileY")

        # Projection of the charge from the deposited to the collected is just multiplying by the CCE
        chgProjProfileX, chgProjProfileY = projectCharge(chgDepProfileX), projectCharge(chgDepProfileY)

        # Diagnostic histograms
        if ((self.logging).level == 10):
            # Create the canvas
            canvas = ROOT.TCanvas("canvas", "Charge collected from the sensors")
            canvas.Divide(2, 2) 
            canvas.cd(1)
            chgDepProfileX.Draw("hist")
            canvas.cd(2)
            chgProjProfileX.Draw("hist")
            canvas.cd(3)
            chgDepProfileY.Draw("hist")
            canvas.cd(4)
            chgProjProfileY.Draw("hist")
            canvas.Update()
            input()
            

        result = {
            'chgDepProfile': (chgDepProfileX.Clone(), chgDepProfileY.Clone()),
            'chgProjProfile': (chgProjProfileX.Clone(), chgProjProfileY.Clone()),
            }
        return result


    # Takes as input data the result of readEdepFromROOT and calculates the uncertainty to attach on each stripfor all the bunches in the run.
    def calculateProfileStatErrs(self, bunchData: dict) -> dict:
        """
        Takes as input data the result of readEdepFromROOT and calculates the uncertainty to attach on each strip for all the bunches in the run.
        
        Parameters
        ----------
            bunchData (dict) : dictionary containing the data of each bunch

        
        Returns
        -------
            npChgDepProfile_err (np.ndarray) : numpy array with the uncertainty on the charge deposited in each strip of each bunch
        
        """
        
        
        # Get the number of bunches
        bunchesNb = len(bunchData)
        
        # Store the profiles in a numpy array
        npChgDepProfile = np.zeros((bunchesNb, 2, 2, 200))         # bunch, det, direction, strip
        
        for bunch in tqdm(bunchData, desc="Calculating errors"):
            # Get the bunch edep 2D Map
            _edepMapUp, _edepMapDo = bunchData[bunch][f"b{bunch}_edepMapUp"], bunchData[bunch][f"b{bunch}_edepMapDo"]

            # Get only the charge deposited profile (0 entry of the result)
            chgDepProfileX_up, chgDepProfileY_up = self.getChgDepProjProfiles(_edepMapUp)['chgDepProfile']
            chgDepProfileX_do, chgDepProfileY_do = self.getChgDepProjProfiles(_edepMapDo)['chgDepProfile']
            
            binsNb = chgDepProfileX_up.GetNbinsX()
            # and store the strip charge into the  numpy array 
            for i in range(binsNb):
                npChgDepProfile[bunch, 0, 0, i] = chgDepProfileX_up.GetBinContent(i);       npChgDepProfile[bunch, 0, 1, i] = chgDepProfileY_up.GetBinContent(i)
                npChgDepProfile[bunch, 1, 0, i] = chgDepProfileX_do.GetBinContent(i);       npChgDepProfile[bunch, 1, 1, i] = chgDepProfileY_do.GetBinContent(i)
        
        # Calculates the strip std between the different bunches
        npChgDepProfile_err = np.zeros((2, 2, 200))     # det, direction, strip
        for det in [0,1]:
            for direction in [0,1]:
                for strip in range(200):
                    stripChgs = npChgDepProfile[:, det, direction, strip]
                    npChgDepProfile_err[det, direction, strip] = np.std(stripChgs)#/np.sqrt(len(stripChgs))
                    

        # Return the errors (useful?)
        return npChgDepProfile_err




    # Read an input ROOT file from the MC simulation and returns the maps of energy deposited in the upstream/downstream sensors
    @dispatch(int, int)
    def readEdepFromROOT(self, bunchParNb: int, pdg: int) -> list:
        """
        Read an input ROOT file from the MC simulation, where the input file is supposed to contain a certain number N*M of physical bunch simulations.
        The function then returns a tuple, with couples of maps of energy deposited in the upstream/downstream sensors.
        Errors are calculated by taking the full statistics and evaluating the std of the strip charge over the population.
        
        Parameters
        ----------
            bunchParNb (int) : number of particles (electrons or gamma) in a bunch
            pdg (int) : PDG code of the particle to consider (default: 22 = gamma) 
            
        Returns
        -------
            list indexed by the bunch nb with per item a dict: { 'bunchParNb': (int) Number particles (electrons or gamma) in the bunch,    "b0_edepMapUp": bunch 0 edepMapUp,  "b0_edepMapDo": bunch 0 edepMapDo }
        """    
        
        # If the _tmp_readEdepFromROOT.npy already exist, then read the data from the file and return immediately
        dumpFname = (self.rFname).replace('.root', f'_bunchParNb{bunchParNb}_pdg{pdg}.npy')
        try:
            msg = f"Loading {dumpFname[dumpFname.rfind('/')+1:]} from file"
            result = np.load(dumpFname, allow_pickle=True)
            msg += " with "+"\x1b[32;1m"+"success"+"\x1b[0m"
            (self.logging).info(msg)
            return result
        except Exception as e:
            msg += f" failed ({e})"
            (self.logging).error(msg)
            (self.logging).info("Calculating from scratch")
        
        # Open the input ROOT file
        riFile = ROOT.TFile.Open(self.rFname, "READ")
        
        # Calculate how many bunches are within this file
        PRIMARY = riFile.ntuple.PRIMARY
        # (pdg == 11) selects only electrons, in future development one may want to consider only positrons or gamma, in which case one has to use (pdg == -11 || pdg == 11) or (pdg == 22)
        primaryParNb = PRIMARY.GetEntries(f"pdg == {pdg}")
        # Number of bunches in the file
        # (bunchChg_pC*1e-12) / (1.160e-19)
        bunchNb = int(primaryParNb / bunchParNb)
        prtTypeStr = f"pdg {pdg}"
        if pdg==11:
            prtTypeStr = "electrons"
        elif pdg==22:
            prtTypeStr = "photons"
        (self.logging).info(f"The number of bunches in the simulation is {bunchNb}. Primary {prtTypeStr} are: {prtTypeStr}")
        
        # Assert that we have non zero bunch number in the sim file
        if bunchNb == 0: raise Exception(f"There are {primaryParNb} primaries and the user number of particles in the bunch is {bunchParNb}. This gets zero particles/bunch with present MC file.")
        if bunchNb > 100: (self.logging).warning(f"Bunch number is large ({bunchNb}). Are you sure this is right?")
        
        # Alias the TTree DUTs
        DUTs = riFile.ntuple.DUTs
        
        # Create the 2D histograms
        edepMapUp_model = ROOT.TH2D("edepMapUp", "Energy deposition in the upstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10)
        edepMapUp_model.GetXaxis().CenterTitle();     edepMapUp_model.GetYaxis().CenterTitle();     edepMapUp_model.GetZaxis().CenterTitle()
        edepMapDo_model = ROOT.TH2D("edepMapDo", "Energy deposition in the downstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10) 
        edepMapDo_model.GetXaxis().CenterTitle();     edepMapDo_model.GetYaxis().CenterTitle();     edepMapDo_model.GetZaxis().CenterTitle()
        
        # Fill the histograms
        result = []
        for i in tqdm(range(bunchNb), desc="Bunch splitting", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed_s:.1f}s<>{remaining_s:.1f}s]'):
            bunch_edepMapUp = edepMapUp_model.Clone(); bunch_edepMapUp.SetName(f"b{i}_edepMapUp")
            bunch_edepMapDo = edepMapDo_model.Clone(); bunch_edepMapDo.SetName(f"b{i}_edepMapDo")
            
            DUTs.Project(f"b{i}_edepMapUp", "edepPosY:edepPosX", f"(edep) * (detID == 0) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
            DUTs.Project(f"b{i}_edepMapDo", "edepPosY:edepPosX", f"(edep) * (detID == 1) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
            
            bunch_edepMapUp.SetDirectory(0); bunch_edepMapDo.SetDirectory(0)        # This way you got the object returned correctly
            
            result.append({ f"b{i}_edepMapUp": bunch_edepMapUp, f"b{i}_edepMapDo": bunch_edepMapDo })
        
        
        # Dump the npChgDepProfile to file for future recalling
        result_np = np.array(result)
        np.save(dumpFname, result_np, allow_pickle=True)
        (self.logging).info(f"[readEdepFromROOT] File {dumpFname[dumpFname.rfind('/')+1:]} saved on file")

        return result_np


    # Overload of readEdepFromROOT to default the pdg to gamma
    @dispatch(int)
    def readEdepFromROOT(self, bunchParNb: int) -> list:
        """
        Read an input ROOT file from the MC simulation, where the input file is supposed to contain a certain number N*M of physical bunch simulations.
        The function then returns a tuple, with couples of maps of energy deposited in the upstream/downstream sensors
        
        Parameters
        ----------
            bunchParNb (int) : number of particles (any type) in a bunch
            
        Returns
        -------
            list where bunch profile is indexed by list index and item contains: { "b0_edepMapUp": bunch 0 edepMapUp,  "b0_edepMapDo": bunch 0 edepMapDo }
        """    
        
        # If the _tmp_readEdepFromROOT.npy already exist, then read the data from the file and return immediately
        bunchDumpFname = (self.rFname).replace('.root', f'_{bunchParNb}.npy')
        if os.path.exists(bunchDumpFname):
            (self.logging).info(f"Loading {bunchDumpFname} from file")
            try:
                # produces something like dummyRun_100k.root -> dummyRun_100k_10_22.npy if 10 is bunchNb and 22 is pdg   
                result = np.load(bunchDumpFname, allow_pickle=True)
                return result[0]
            except Exception as e:
                (self.logging).error(f"Loading failed due to {e}")
        
        
        # Open the input ROOT file
        riFile = ROOT.TFile.Open(self.rFname, "READ")
        
        # Calculate how many bunches are within this file
        PRIMARY = riFile.ntuple.PRIMARY
        # (pdg == 11) selects only electrons, in future development one may want to consider only positrons or gamma, in which case one has to use (pdg == -11 || pdg == 11) or (pdg == 22)
        primaryParNb = PRIMARY.GetEntries()
        # Number of bunches in the file
        # (bunchChg_pC*1e-12) / (1.160e-19)
        bunchNb = int(primaryParNb / bunchParNb)
        (self.logging).debug(f"The number of bunches in the simulation is {bunchNb}. Primary are: {primaryParNb}")
        
        # Assert that we have non zero bunch number in the sim file
        if bunchNb == 0: raise Exception(f"There are {primaryParNb} primaries and the user number of particles in the bunch is {bunchParNb}. This gets zero particles/bunch with present MC file.")
        # Send warning if the number of bunches is large
        if bunchNb > 100: (self.logging).warning(f"Bunch number is large ({bunchNb}). Are you sure this is right?")
        
        # Alias the TTree DUTs
        DUTs = riFile.ntuple.DUTs
        # Create the 2D histograms
        edepMapUp_model = ROOT.TH2D("edepMapUp", "Energy deposition in the upstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10)
        edepMapUp_model.GetXaxis().CenterTitle();     edepMapUp_model.GetYaxis().CenterTitle();     edepMapUp_model.GetZaxis().CenterTitle()
        edepMapDo_model = ROOT.TH2D("edepMapDo", "Energy deposition in the downstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10) 
        edepMapDo_model.GetXaxis().CenterTitle();     edepMapDo_model.GetYaxis().CenterTitle();     edepMapDo_model.GetZaxis().CenterTitle()
        
        # Fill the histograms
        result = []
        for i in tqdm(range(bunchNb), desc="Bunch splitting"):
            bunch_edepMapUp = edepMapUp_model.Clone(); bunch_edepMapUp.SetName(f"b{i}_edepMapUp")
            bunch_edepMapDo = edepMapDo_model.Clone(); bunch_edepMapDo.SetName(f"b{i}_edepMapDo")
            
            DUTs.Project(f"b{i}_edepMapUp", "edepPosY:edepPosX", f"(edep) * (detID == 0) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
            DUTs.Project(f"b{i}_edepMapDo", "edepPosY:edepPosX", f"(edep) * (detID == 1) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
            
            bunch_edepMapUp.SetDirectory(0); bunch_edepMapDo.SetDirectory(0)        # This way you got the object returned correctly
            
            result.append({ f"b{i}_edepMapUp": bunch_edepMapUp, f"b{i}_edepMapDo": bunch_edepMapDo })
        
        
        # Dump the npChgDepProfile to file for future recalling
        result_np = np.array(result)
        np.save(bunchDumpFname, result_np, allow_pickle=True)
        (self.logging).info(f"Dumped {bunchDumpFname} to file")
        
        return result_np


    # Read an input ROOT file from the MC simulation and returns the maps of energy deposited in the upstream/downstream sensors
    @dispatch()
    def readEdepFromROOT(self) -> list:
        return self.readEdepFromROOT(self.rFname, 1)


    #################################################################################################

    @dispatch(int, float, float)
    def readProjProfilesFromROOT(self, bunchParNb: int, cce: float, avgPairEn: float) -> np.array:
        """
        Loads the MC data from the fname, slice the total energy deposition in the sensor into bunches with 'bunchParNb' particles each,
        and returns a collection of TGraphErrors with the projected charge profiles (cce and avgPairEn are used for this) whose strip charges have also
        errors attached to them. The errors are calculated in the following way:
        1. the strip charge i is evaluated for each bunch number and the std is calculated,
        2. to the strip charge i of each bunch is then attached this std value.
        
        Parameters
        ----------
            bunchParNb (int) : number of particles (electrons or gamma) in a bunch
            cce (float) : Charge collection efficiency of the sensor for the bare geometrical projection of the dep chg. to proj chg.
            avgPairEn (float) : Average energy required to create an electron-hole pair in the sensor (in eV)
        
        Returns
        -------
            writeme (np.array) : text   
        """
        
        # Attempt to load from file
        dumpFname = (self.rFname).replace('.root', f'_bunchParNb{bunchParNb}_cce{cce}_avgPairEn{avgPairEn}.npy')    
        if os.path.exists(dumpFname):
            msg = f"Loading {dumpFname[dumpFname.rfind('/')+1:]} from file"
            try:
                depProjProfsBunch_np = np.load(dumpFname, allow_pickle=True)
                msg += " with "+"\x1b[32;1m"+"success"+"\x1b[0m"
                (self.logging).info(msg)
                
                # Debug
                if ((self.logging).level == 10):
                    depProjProfsBunch_np[0]['d0_x'].Draw("AP")
                    input()
                return depProjProfsBunch_np
            except Exception as e:
                msg += f" failed ({e})"
                (self.logging).error(msg)
                (self.logging).info("Calculating from scratch")
        

        # Get - from the ROOT file produced by the MC - the list with, for each bunch, the map of energy deposition
        bunchEMaps = self.readEdepFromROOT(bunchParNb)
        # bunchEMaps is a list like [obj1, obj2, ..., objN ] with N the number of bunches and objI made like
        # this objI = { f"b{i}_edepMapUp": bunch_edepMapUp, f"b{i}_edepMapDo": bunch_edepMapDo }
        
        # Calculate the statistical error to be attached to each strip
        npChgDepProfile_err = self.calculateProfileStatErrs(bunchEMaps)
        
        # Attach the errors to the profiles
        depProjProfsBunch = []
        for bunch in tqdm(bunchEMaps, desc="ProjGraphs with err"):
            # Get the bunch edep 2D Map
            _edepMapUp, _edepMapDo = bunchEMaps[bunch][f"b{bunch}_edepMapUp"], bunchEMaps[bunch][f"b{bunch}_edepMapDo"]
            
            # Get the charge projected profile upstream downstream in both directions
            chgProjProfileX_up, chgProjProfileY_up = self.getChgDepProjProfiles(_edepMapUp)['chgProjProfile']
            chgProjProfileX_do, chgProjProfileY_do = self.getChgDepProjProfiles(_edepMapDo)['chgProjProfile']
            
            # Attach the errors to the bunch profiles
            for i in range(200):
                chgProjProfileX_up.SetBinError(i, (npChgDepProfile_err[0, 0, i] * cce) )
                #chgProjProfileX_do.SetBinError(i, (npChgDepProfile_err[0, 1, i] * cce) )
                #chgProjProfileY_up.SetBinError(i, (npChgDepProfile_err[1, 0, i] * cce) )
                chgProjProfileY_do.SetBinError(i, (npChgDepProfile_err[1, 1, i] * cce) )
            
            # Create a TGraph with errors and attach the profile of projected charge from bunch i
            depProjProfsBunch.append({ 'd0_x' : ROOT.TGraphErrors(chgProjProfileX_up), 'd1_y' : ROOT.TGraphErrors(chgProjProfileY_do) })
        

        # Remove the objects from memory
        (self.logging).critical("Remove implement me")
        #for tp in ['Dep', 'Proj']:
        #    for dir in ['X', 'Y']:
        #        obName = f'chg{tp}Profile{dir}'
        #        ptr = ROOT.Get(obName)
        #        ptr.Delete()
        #        del ptr
            
        
        # Dump the processed profiles in memory into a file
        depProjProfsBunch_np = np.array(depProjProfsBunch)
        np.save(dumpFname, depProjProfsBunch_np, allow_pickle=True)
        (self.logging).info(f"[readProjProfilesFromROOT] File {dumpFname[dumpFname.rfind('/')+1:]} saved on file")
        
        
        if ((self.logging).level == 10):
            # Debug
            depProjProfsBunch_np[0]['d0_x'].Draw("AP")
            input()
        
        return depProjProfsBunch_np
    
    
    @dispatch()
    def readProjProfilesFromROOT(self) -> np.array:
        return self.readProjProfilesFromROOT(self.bunchParNb, self.cce, self.avgPairEn)



    @dispatch(int, int, int, float, float)
    def readProjProfilesFromROOT(self, bunchParNb: int, pdg: int, cce: float, avgPairEn: float) -> np.array:
        (self.logging).error("I'm not implemented yet")
        raise Exception("I'm not implemented yet")
    
    
    
    

# Test function for the class
if __name__ == "__main__":
    # Open the dummyFile and plot profile of the first bunch
    rMCClass = readFromMc("build/dummyRun_100k.root", 10000, 0.2, 27.0)
    # Enable debug
    #rMCClass.logging.setLevel(10)
    print(rMCClass.readProjProfilesFromROOT())
    exit()