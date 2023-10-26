# Description: This script contains the functions to project and digitize the energy deposition in the sensors.
#              The projection is done by projecting the energy deposition in the sensors to the strips.
#              The digitization is done in 3 steps:
#              1. Apply the front-end noise to the profile of charge collected at the strips
#              2. Apply the front-end amplification to the profile of charge projected at the strips
#              3. Convert a profile from continuous voltage value to discrete value
from multipledispatch import dispatch
from matplotlib import pyplot as plt
from tqdm import tqdm
import numpy as np
import ROOT
from rdataStruct import rdataStruct_OPT


# digitization.py logger
from logger import create_logger
logging = create_logger("digitization")
    


# Enable ROOT implicit multithreading
ROOT.EnableImplicitMT(16)


# Note:
# The number of entries in a weighted histogram are the rounded sum of the weights! Be aware! 

# Parameters     
_stripPitch = 0.1  
_nbStrips = 192   

# Physical constants 
_cce = 0.2
_avgPairCreationEnergy = 27.0 # in eV


# Read an input ROOT file from the MC simulation and returns the maps of energy deposited in the upstream/downstream sensors
@dispatch(str)
def readEdepFromROOT(fname: str) -> tuple:
    """
    Read an input ROOT file from the MC simulation and returns the maps of energy deposited in the upstream/downstream sensors
    
    Parameters
    ----------
        fname (str) : filename of the ROOT file produced from the StandaloneMc simulation
        
    Returns
    -------
        edepMapUp (ROOT.TH2D) : 2D ROOT histogram with energy deposition in the upstream detector 
        edepMapDo (ROOT.TH2D) : 2D ROOT histogram with energy deposition in the downstream detector 
    """
    
    # Open the input ROOT file
    riFile = ROOT.TFile.Open(fname, "READ")
    # Alias the TTree DUTs
    DUTs = riFile.ntuple.DUTs
    
    # Create the 2D histograms
    edepMapUp = ROOT.TH2D("edepMapUp", "Energy deposition in the upstream sensor; X [mm]; Y [mm]; energy [keV]", 200, -10, 10, 200, -10, 10)
    edepMapUp.GetXaxis().CenterTitle()
    edepMapUp.GetYaxis().CenterTitle()
    edepMapUp.GetZaxis().SetTitle("Energy deposited [keV]")
    edepMapUp.GetZaxis().CenterTitle()

    edepMapDo = ROOT.TH2D("edepMapDo", "Energy deposition in the downstream sensor; X [mm]; Y [mm]; energy [keV]", 200, -10, 10, 200, -10, 10) 
    edepMapDo.GetXaxis().CenterTitle()
    edepMapDo.GetYaxis().CenterTitle()
    edepMapDo.GetZaxis().SetTitle("Energy deposited [keV]")
    edepMapDo.GetZaxis().CenterTitle()

    # Fill the histograms
    DUTs.Project("edepMapUp", "edepPosY:edepPosX", "edep * (detID == 0)")
    DUTs.Project("edepMapDo", "edepPosY:edepPosX", "edep * (detID == 1)")
    
    ## Debug
    #canvas = ROOT.TCanvas("canvas", "canvas")
    #canvas.Divide(2,1)
    #canvas.cd(1)
    #edepMapUp.Draw("colz")
    #canvas.cd(2)
    #edepMapDo.Draw("colz")
    #input()
    #exit()
    

    # This way you got the object returned correctly
    edepMapUp.SetDirectory(0)
    edepMapDo.SetDirectory(0)
    return (edepMapUp, edepMapDo)


# Read an input ROOT file from the MC simulation and returns the maps of energy deposited in the upstream/downstream sensors
@dispatch(str, int, int)
def readEdepFromROOT(fname: str, bunchParNb: int, pdg: int) -> dict:
    """
    Read an input ROOT file from the MC simulation, where the input file is supposed to contain a certain number N*M of physical bunch simulations.
    The function then returns a tuple, with couples of maps of energy deposited in the upstream/downstream sensors
    
    Parameters
    ----------
        fname (str) : filename of the ROOT file produced from the Standalone MC simulation
        bunchParNb (int) : number of particles (electrons or gamma) in a bunch
        pdg (int) : PDG code of the particle to consider (default: 22 = gamma) 
        
    Returns
    -------
        { bunchID:
            { 'bunchParNb': (int) Number particles (electrons or gamma) in the bunch,    "b0_edepMapUp": bunch 0 edepMapUp,  "b0_edepMapDo": bunch 0 edepMapDo }
        }
    """    
    
    # If the _tmp_readEdepFromROOT.npy already exist, then read the data from the file and return immediately
    bunchDumpFname = fname.replace('.root', f'_{bunchParNb}_{pdg}.npy')
    try:
        # produces something like dummyRun_100k.root -> dummyRun_100k_10_22.npy if 10 is bunchNb and 22 is pdg   
        print(f"Loading {bunchDumpFname} from file")
        result = np.load(bunchDumpFname, allow_pickle=True)
        return result[0]
    except:
        pass
    
    # Open the input ROOT file
    riFile = ROOT.TFile.Open(fname, "READ")
    
    # Calculate how many bunches are within this file
    PRIMARY = riFile.ntuple.PRIMARY
    # (pdg == 11) selects only electrons, in future development one may want to consider only positrons or gamma, in which case one has to use (pdg == -11 || pdg == 11) or (pdg == 22)
    primaryParNb = PRIMARY.GetEntries(f"pdg == {pdg}")
    # Number of bunches in the file
    # (bunchChg_pC*1e-12) / (1.160e-19)
    bunchNb = int(primaryParNb / bunchParNb)
    print(f"The number of bunches in the simulation is {bunchNb}.")
    prtTypeStr = f"pdg {pdg}"
    if pdg==11:
        prtTypeStr = "electrons"
    elif pdg==22:
        prtTypeStr = "photons"
    print("Primary "+prtTypeStr+" are:", primaryParNb)
    
    # Assert that we have non zero bunch number in the sim file
    if bunchNb == 0: raise Exception(f"There are {primaryParNb} primaries and the user number of particles in the bunch is {bunchParNb}. This gets zero particles/bunch with present MC file.")
    
    if bunchNb > 100: print(f"Bunch number is large ({bunchNb}). Are you sure this is right?")
    
    # Alias the TTree DUTs
    DUTs = riFile.ntuple.DUTs
    
    # Create the 2D histograms
    edepMapUp_model = ROOT.TH2D("edepMapUp", "Energy deposition in the upstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10)
    edepMapUp_model.GetXaxis().CenterTitle();     edepMapUp_model.GetYaxis().CenterTitle();     edepMapUp_model.GetZaxis().CenterTitle()
    edepMapDo_model = ROOT.TH2D("edepMapDo", "Energy deposition in the downstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10) 
    edepMapDo_model.GetXaxis().CenterTitle();     edepMapDo_model.GetYaxis().CenterTitle();     edepMapDo_model.GetZaxis().CenterTitle()
    
    # Fill the histograms
    result = {}
    for i in tqdm(range(bunchNb), desc="Bunch splitting", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed_s:.1f}s<>{remaining_s:.1f}s]'):
        bunch_edepMapUp = edepMapUp_model.Clone(); bunch_edepMapUp.SetName(f"b{i}_edepMapUp")
        bunch_edepMapDo = edepMapDo_model.Clone(); bunch_edepMapDo.SetName(f"b{i}_edepMapDo")
        
        DUTs.Project(f"b{i}_edepMapUp", "edepPosY:edepPosX", f"(edep) * (detID == 0) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
        DUTs.Project(f"b{i}_edepMapDo", "edepPosY:edepPosX", f"(edep) * (detID == 1) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
        
        bunch_edepMapUp.SetDirectory(0); bunch_edepMapDo.SetDirectory(0)        # This way you got the object returned correctly
        
        result[i] = { 'bunchNb': bunchNb, f"b{i}_edepMapUp": bunch_edepMapUp, f"b{i}_edepMapDo": bunch_edepMapDo }
    
    
    # Dump the npChgDepProfile to file for future recalling
    np.save(bunchDumpFname, np.array([result]), allow_pickle=True)
    print(f"Dumped {bunchDumpFname} to file")
    
    return result


#  Overload of readEdepFromROOT to default the pdg to gamma
@dispatch(str, int)
def readEdepFromROOT(fname: str, bunchParNb: int) -> dict:
    """
    Read an input ROOT file from the MC simulation, where the input file is supposed to contain a certain number N*M of physical bunch simulations.
    The function then returns a tuple, with couples of maps of energy deposited in the upstream/downstream sensors
    
    Parameters
    ----------
        fname (str) : filename of the ROOT file produced from the Standalone MC simulation
        bunchParNb (int) : number of particles (any type) in a bunch
        
    Returns
    -------
        { bunchID:
            { 'bunchParNb': (int) Number particles (any type) in the bunch,    "b0_edepMapUp": bunch 0 edepMapUp,  "b0_edepMapDo": bunch 0 edepMapDo }
        }
    """    
    
    # If the _tmp_readEdepFromROOT.npy already exist, then read the data from the file and return immediately
    bunchDumpFname = fname.replace('.root', f'_{bunchParNb}.npy')
    try:
        # produces something like dummyRun_100k.root -> dummyRun_100k_10_22.npy if 10 is bunchNb and 22 is pdg   
        print(f"Loading {bunchDumpFname} from file")
        result = np.load(bunchDumpFname, allow_pickle=True)
        return result[0]
    except:
        pass
    
    # Open the input ROOT file
    riFile = ROOT.TFile.Open(fname, "READ")
    
    # Calculate how many bunches are within this file
    PRIMARY = riFile.ntuple.PRIMARY
    # (pdg == 11) selects only electrons, in future development one may want to consider only positrons or gamma, in which case one has to use (pdg == -11 || pdg == 11) or (pdg == 22)
    primaryParNb = PRIMARY.GetEntries()
    # Number of bunches in the file
    # (bunchChg_pC*1e-12) / (1.160e-19)
    bunchNb = int(primaryParNb / bunchParNb)
    print(f"The number of bunches in the simulation is {bunchNb}.")
    print("Primary are:", primaryParNb)
    
    # Assert that we have non zero bunch number in the sim file
    if bunchNb == 0: raise Exception(f"There are {primaryParNb} primaries and the user number of particles in the bunch is {bunchParNb}. This gets zero particles/bunch with present MC file.")
    
    if bunchNb > 100: print(f"Bunch number is large ({bunchNb}). Are you sure this is right?")
    
    # Alias the TTree DUTs
    DUTs = riFile.ntuple.DUTs
    
    # Create the 2D histograms
    edepMapUp_model = ROOT.TH2D("edepMapUp", "Energy deposition in the upstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10)
    edepMapUp_model.GetXaxis().CenterTitle();     edepMapUp_model.GetYaxis().CenterTitle();     edepMapUp_model.GetZaxis().CenterTitle()
    edepMapDo_model = ROOT.TH2D("edepMapDo", "Energy deposition in the downstream sensor; X [mm]; Y [mm]; Energy deposited [keV]", 200, -10, 10, 200, -10, 10) 
    edepMapDo_model.GetXaxis().CenterTitle();     edepMapDo_model.GetYaxis().CenterTitle();     edepMapDo_model.GetZaxis().CenterTitle()
    
    # Fill the histograms
    result = {}
    for i in tqdm(range(bunchNb), desc="Bunch splitting", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed_s:.1f}s<>{remaining_s:.1f}s]'):
        bunch_edepMapUp = edepMapUp_model.Clone(); bunch_edepMapUp.SetName(f"b{i}_edepMapUp")
        bunch_edepMapDo = edepMapDo_model.Clone(); bunch_edepMapDo.SetName(f"b{i}_edepMapDo")
        
        DUTs.Project(f"b{i}_edepMapUp", "edepPosY:edepPosX", f"(edep) * (detID == 0) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
        DUTs.Project(f"b{i}_edepMapDo", "edepPosY:edepPosX", f"(edep) * (detID == 1) * (event >= {bunchParNb * i} && event < {bunchParNb * (i+1)})")
        
        bunch_edepMapUp.SetDirectory(0); bunch_edepMapDo.SetDirectory(0)        # This way you got the object returned correctly
        
        result[i] = { 'bunchNb': bunchNb, f"b{i}_edepMapUp": bunch_edepMapUp, f"b{i}_edepMapDo": bunch_edepMapDo }
    
    
    # Dump the npChgDepProfile to file for future recalling
    np.save(bunchDumpFname, np.array([result]), allow_pickle=True)
    print(f"Dumped {bunchDumpFname} to file")
    
    return result





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
    chgProjProfile.Scale(_cce)
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

# Takes the map of energy deposited in a sensor and returns the map of charge collected at the strips (hor./vert.)
def getChgProfiles(edepHist: ROOT.TH2D) -> dict:
    """
    Takes the map of energy deposited in a sensor and returns the map of charge collected at the strips (hor./vert.)
    
    Parameters
    ----------
        edepHist (ROOT.TH2D) : histogram with the map of energy deposition in a sensor
        
    Returns
    -------
        'chgDepProfile': (chgDepProfileX.Clone(), chgDepProfileY.Clone()),
        'chgProjProfile': (chgProjProfileX.Clone(), chgProjProfileY.Clone())
    """

    # Get horizontal and vertical profiles of energy depositions in the sensor
    edepProfileX, edepProfileY = edepHist.ProjectionX(), edepHist.ProjectionY()
        
    # This maps the profile of edeps to a profile of charge deposited (in C)
    # The 1e3 accounts for the energy in the histo being in keV
    eCharge = 1.602176634e-19
    chgDepProfileX, chgDepProfileY = edepProfileX.Clone(), edepProfileY.Clone()
    chgDepProfileX.Scale(1e3 * eCharge /_avgPairCreationEnergy), chgDepProfileY.Scale(1e3 * eCharge /_avgPairCreationEnergy)

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
    if (logging.level == 10):
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
def computeError_v2(bunchData: dict) -> dict:
    """
    Takes as input data the result of readEdepFromROOT and calculates the uncertainty to attach on each stripfor all the bunches in the run.
    
    Parameters
    ----------
    bunchData (dict) : dictionary containing the data of each bunch
    
    Returns
    -------
    npChgDepProfile_err (np.ndarray) : numpy array with the uncertainty on the charge deposited in each strip of each bunch
    
    """
    
    # If the _tmp_npChgDepProfile_err.npy already exist, then read the data from the file and return immediately
    try:
        print("loading npChgDepProfile_err from file")              # GENERALIZE TO MATCH THE BUNCH DATA FILE
        npChgDepProfile_err = np.load("_tmp_npChgDepProfile_err.npy", allow_pickle=True)
        #plt.plot(npChgDepProfile_err[0,0,:], marker='o', linestyle='', )
        #plt.show()
        return npChgDepProfile_err
    except:
        pass
    
    
    # Get the number of bunches
    bunchesNb = len(bunchData)
    
    # Store the profiles in a numpy array
    npChgDepProfile = np.zeros((bunchesNb, 2, 2, 200))         # bunch, det, direction, strip
    
    for bunch in bunchData:
        # Get the bunch edep 2D Map
        _edepMapUp, _edepMapDo = bunchData[bunch][f"b{bunch}_edepMapUp"], bunchData[bunch][f"b{bunch}_edepMapDo"]

        # Get only the charge deposited profile (0 entry of the result)
        chgDepProfileX_up, chgDepProfileY_up = getChgProfiles(_edepMapUp)['chgDepProfile']
        chgDepProfileX_do, chgDepProfileY_do = getChgProfiles(_edepMapDo)['chgDepProfile']
        
        binsNb = chgDepProfileX_up.GetNbinsX()
        # and store the strip charge into the  numpy array 
        for i in range(binsNb):
            npChgDepProfile[bunch, 0, 0, i] = chgDepProfileX_up.GetBinContent(i);       npChgDepProfile[bunch, 0, 1, i] = chgDepProfileY_up.GetBinContent(i)
            npChgDepProfile[bunch, 1, 0, i] = chgDepProfileX_do.GetBinContent(i);       npChgDepProfile[bunch, 1, 1, i] = chgDepProfileY_do.GetBinContent(i)
    
    
    npChgDepProfile_err = np.zeros((2, 2, 200))
    for det in [0,1]:
        for direction in [0,1]:
            for strip in range(200):
                stripChgs = npChgDepProfile[:, det, direction, strip]
                npChgDepProfile_err[det, direction, strip] = np.std(stripChgs)/np.sqrt(len(stripChgs))
    
    #
    # Dump the npChgDepProfile to file for future recalling
    np.save("_tmp_npChgDepProfile_err.npy", npChgDepProfile_err, allow_pickle=True)
    
    #
    plt.errorbar(np.linspace(0,200,200), npChgDepProfile[0,0,0,:], xerr=0, yerr=npChgDepProfile_err[0,0,:], marker='o', linestyle='', capsize=5, label='Data with Error Bars')
    plt.show()
    return npChgDepProfile




# From the 2D maps of energy depositions in the upstream/downstream sensors it produces a plot with
# 1. the profiles of charge deposited in the sensors
# 2. the profiles of charge projected at the strips after applying the projectCharge function
def processPlot_chgDepProj(edepMapUp: ROOT.TH2D, edepMapDo: ROOT.TH2D) -> None:
    """
    From the 2D maps of energy depositions in the upstream/downstream sensors it produces a plot with
    1. the profiles of charge deposited in the sensors
    2. the profiles of charge projected at the strips after applying the projectCharge function

    Parameters
    ----------
        edepMapUp (ROOT.TH2D) : histogram with the map of energy deposition in the upstream sensor
        edepMapDo (ROOT.TH2D) : histogram with the map of energy deposition in the downstream sensor    
    """
    
    # Get the profiles of charge deposited and projected
    chgProfUp, chgProfDo = getChgProfiles(edepMapUp), getChgProfiles(edepMapDo)

    # Create the canvas
    canvas = ROOT.TCanvas("canvas", "Charge collected from the sensors")
    canvas.Divide(2, 2) 
    canvas.cd(1)
    chgProfUp['chgDepProfile'][0].Draw("hist")
    canvas.cd(2)
    chgProfUp['chgProjProfile'][0].Draw("hist")
    canvas.cd(3)
    chgProfDo['chgDepProfile'][1].Draw("hist")
    canvas.cd(4)
    chgProfDo['chgProjProfile'][1].Draw("hist")
    canvas.Update()
    
    # Save the canvas to pdf for showing to M. Morandin
    canvas.SaveAs("chgDepProj.pdf")


# Converts the X/Y position to a strip number between 1-192. By convention, strip #1 is the one with negative X/Y coordinate
def calculateStripNb(pos_mm: float) -> int:
    """
    Converts the X/Y position to a strip number between 1-192. By convention, strip #1 is the one with negative X/Y coordinate
    
    Parameters
    ----------
        pos_mm (float) : position in mm
        
    Returns
    -------
        stripNb (int) : strip number between 1-192
    """
    
    return int( (pos_mm + _stripPitch*_nbStrips/2.0) / (_stripPitch*_nbStrips) * _nbStrips + 1.0 )
    


# Change labels to have strip numbers in the horizontal axis (instead of X/Y in mm)
def mapStripProfile(hist: ROOT.TH1D) -> ROOT.TH1D:
    """
    Change labels to have strip numbers in the horizontal axis (instead of X/Y in mm)
    
    Parameters
    ----------
        hist (ROOT.TH1D) : input histogram
        
    Returns
    -------
        newHist (ROOT.TH1D) : histogram with strip numbers in the horizontal axis
    """
    
    newHist = ROOT.TH1D(hist.GetName()+'_strip', hist.GetTitle(), _nbStrips, 1, _nbStrips)
    newHist.GetXaxis().SetTitle("strip number [1-192]")
    newHist.GetXaxis().CenterTitle()
    newHist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    newHist.GetYaxis().CenterTitle()
    # Check if the histogram is 2D
    if hist.GetDimension() == 2:
    # If so, get a 1D projection along the x-axis
        hist = hist.ProjectionX()

    for bin in range(1, hist.GetNbinsX()+1):
        strip = calculateStripNb(hist.GetBinCenter(bin))
        bin_content = hist.GetBinContent(bin)
        newHist.Fill(strip, bin_content)
        #bin_idx = newHist.FindBin(strip)
        #newHist.SetBinError(bin_idx, dstrip)
    return newHist



# Simulates the cross talk between two strips of the sensors. The cross talk is simulated by adding a normal distribution centered around the selected strip number for each sensor.
@dispatch(ROOT.TH1D, list)
def simulateCrosstalk(chgProfile: ROOT.TH1D, chgShrNN: list) -> ROOT.TH1D:
    """
    Simulates the cross-talk between neighbour strips. The cross talk is simulated by adding a fraction of the charge of each strip to the nearest neighbours, with
    a perncentage depending on the diostance between strip i and strip j

    Parameters
    ----------
        chgProfile (ROOT.TH1D) : input histogram with the charge profile of the sensor
        chgShrNN (list) : list of fractions of charge to share with the neighbours
        
    Returns
    -------
        chgProfileWithCrossTalk (ROOT.TH1D) : histogram with the charge profile of the sensor after the cross-talk simulation
    """
    stripNb = chgProfile.GetNbinsX()
    npProfile = np.zeros((stripNb,2))
    for i in range(1, stripNb+1):
        # Get the charge and error on it
        stripChg = chgProfile.GetBinContent(i)
        stripChg_err = chgProfile.GetBinError(i)
        npProfile[i-1, 0] = stripChg
        npProfile[i-1, 1] = stripChg_err

    # Create arrays to store the crosstalk contributions and their errors
    crosstalk_contribs = np.zeros((stripNb, len(chgShrNN)))
    crosstalk_contribs_err = np.zeros((stripNb, len(chgShrNN)))

    # Calculate the crosstalk contributions
    for i in range(stripNb):
        if npProfile[i, 0] != 0:
            for j in range(len(chgShrNN)):
                crosstalk_contribs[i, j] = npProfile[i, 0] * chgShrNN[j]
                crosstalk_contribs_err[i, j] = npProfile[i, 1] * chgShrNN[j]

    # Update the charge profile
    updatedProfile = npProfile.copy()
    #initialChg = np.sum(updatedProfile, 0)
    for i in range(stripNb):
        if npProfile[i, 0] != 0:
            for j in range(len(chgShrNN)):
                # Update the charge on the left neighbor
                if i-(j+1) >= 0:
                    updatedProfile[i-(j+1), 0] += crosstalk_contribs[i, j]
                    updatedProfile[i-(j+1), 1] = np.sqrt(updatedProfile[i-(j+1), 1]**2 + crosstalk_contribs_err[i, j]**2)
                
                # Update the charge on the right neighbor
                if i+(j+1) < stripNb:
                    updatedProfile[i+(j+1), 0] += crosstalk_contribs[i, j]
                    updatedProfile[i+(j+1), 1] = np.sqrt(updatedProfile[i+(j+1), 1]**2 + crosstalk_contribs_err[i, j]**2)
                
                # Update the charge on current strip  
                updatedProfile[i, 0] -= 2*crosstalk_contribs[i, j]
                updatedProfile[i, 1] = 2*np.sqrt(updatedProfile[i, 1]**2 + crosstalk_contribs_err[i, j]**2)
                
    ## Normalize the total charge shared
    #currentChg = np.sum(updatedProfile, 0)
    #updatedProfile *=  initialChg/currentChg
    
    # Diagnostic plots
    if (logging.level == 10):
        view, ax = plt.subplots()
        view.suptitle("Projected charge profile before/after crosstalk charge spread")
        ax.set_xlabel("strip no.")
        ax.set_ylabel("strip charge [C]")

        # Use the same 'x' data for both "before" and "after" plots
        x_data = np.arange(1, stripNb + 1)

        ax.errorbar(x_data, npProfile[:, 0], yerr=npProfile[:, 1], fmt='o', label='before crosstalk')
        ax.errorbar(x_data, updatedProfile[:, 0], yerr=updatedProfile[:, 1], fmt='o', label='after crosstalk')
        ax.legend(loc="upper right")
        view.show()
        plt.show()


    # Create a new histogram for the updated charge profile
    chgProfileWithCrossTalk = chgProfile.Clone()
    chgProfileWithCrossTalk.SetName(chgProfileWithCrossTalk.GetName().replace("Proj", "ProjWithCrosstalk"))
    chgProfileWithCrossTalk.SetTitle(chgProfileWithCrossTalk.GetTitle() + " with cross-talk applied")
    chgProfileWithCrossTalk.GetYaxis().SetTitle(chgProfileWithCrossTalk.GetYaxis().GetTitle().replace("projected", "projected with CT"))
    
    for i in range(1, stripNb+1):
        chgProfileWithCrossTalk.SetBinContent(i, updatedProfile[i-1, 0])
        #chgProfileWithCrossTalk.SetBinError(i, updatedProfile[i-1, 1])
    return chgProfileWithCrossTalk


# Overload to accomodate the structure of the other functions
@dispatch(dict, dict, list)
def simulateCrosstalk(chgProfilesUp: dict, chgProfilesDo: dict, chgShrNN: list) -> tuple:
    """
    Returns the tuple with the dictionaries chgProfilesUp and chgProfilesDo with the new 'chgProjWithCTProfile' key
    containing the profiles with the strip cross-talk effect applied
    
    Parameters
    ----------
        chgProfilesUp (dict) : dictionary with the charge profiles of the upstream sensor
        chgProfilesDo (dict) : dictionary with the charge profiles of the downstream sensor
        chgShrNN (list) : list of fractions of charge to share with the neighbours
        
    Returns
    -------
        (chgProfilesUp, chgProfilesDo) (tuple) : tuple with the dictionaries chgProfilesUp and chgProfilesDo with the new 'chgProjWithCTProfile' key
    """
    chgProjWithCTProfileUp_X = simulateCrosstalk(chgProfilesUp['chgProjProfile'][0], chgShrNN)
    chgProjWithCTProfileUp_Y = simulateCrosstalk(chgProfilesUp['chgProjProfile'][1], chgShrNN)
    chgProfilesUp['chgProjWithCTProfile'] = (chgProjWithCTProfileUp_X, chgProjWithCTProfileUp_Y)
    
    chgProjWithCTProfileDo_X = simulateCrosstalk(chgProfilesUp['chgProjProfile'][0], chgShrNN)
    chgProjWithCTProfileDo_Y = simulateCrosstalk(chgProfilesUp['chgProjProfile'][1], chgShrNN)
    chgProfilesDo['chgProjWithCTProfile'] = (chgProjWithCTProfileDo_X, chgProjWithCTProfileDo_Y)
    
    return (chgProfilesUp, chgProfilesDo)
    
    
    
    
    

################################################################################################################################################

# Take an 2D histogram of energy deposited and returns a numpy array with bin position and energy deposited (x,y,edep)
def convertROOT2DHist_toNumpy(edepHist : ROOT.TH2D) -> np.array:
    """
    Take an 2D histogram of energy deposited and returns a numpy array with bin position and energy deposited (x,y,edep)
    
    Parameters
    ----------
        hist (ROOT.TH2D) : input histogram with the 2D map of energy deposition
        
    Returns
    -------
        np.array : numpy array with bin position and energy deposited (x,y,edep)
    """
    
    binx = edepHist.GetNbinsX()
    biny = edepHist.GetNbinsY()
    if binx != biny: raise Exception(f"Unexpected different number of bins in X and Y axes: {binx} | {biny}")
    # From here on you are sure that buins are the same numbers
    
    # Create a multidimensional (3) numpy array with x,y,content
    binXCenters = np.zeros(binx)
    binYCenters = np.zeros(binx)
    binContents = np.zeros(binx)
    
    for i in range(binx):
        binXCenters[i] = edepHist.GetXaxis().GetBinCenter(i)
        binYCenters[i] = edepHist.GetYaxis().GetBinCenter(i)
        binContents[i] = edepHist.GetBinContent(i, i)

    return np.column_stack((binXCenters, binYCenters, binContents))



# Map all histograms in the given dictionary using mapStripProfile.
def mapAllHistograms(histogram_dict: dict) -> dict:
    """
    Map all histograms in the given dictionary using mapStripProfile.
    
    Parameters
    ----------
        histogram_dict (dict) : dictionary of ROOT.TH1D histograms
        
    Returns
    -------
        new_dict (dict) : dictionary of mapped ROOT.TH1D histograms
    """

    new_dict = {}
    for key, hist in histogram_dict.items():
        new_hist = mapStripProfile(hist)

        # Scale all histograms
        new_hist.Scale(1e-14)

        new_dict[key] = new_hist

    return new_dict

################################################################################################################################################

# Dump all the objects (MUST BE ROOT OBJECT, BE AWARE) in the args to the ROOT file with filename 'fname'
def dumpAllInRoot(fname: str, *args):
    """
    Dump all the objects (MUST BE ROOT OBJECT, BE AWARE) in the args to the ROOT file with filename 'fname'
    """
    print("All the objects will be saved into: ", fname)
    roFile = ROOT.TFile(fname, "RECREATE")
    
    for arg in args:
        # Check if the arg is a ROOT histogram
        if isinstance(arg, (ROOT.TH1D, ROOT.TH2D)):
            # If it is, write it to the file
            print(f"Dumping histogram {arg.GetName()} ...", end='')
            arg.Write(f"hist_{arg.GetName()}")
            print("done")
        # If the arg is a dictionary
        elif isinstance(arg, dict): 
            # Iterate over its items
            for key, hist in arg.items():
                print(f"Dumping histogram {key} ...", end='')
                hist.Write(f"hist_{key}")
                print("done")
        else:
            print(f"Unsupported arg type: {type(arg)}")

    roFile.Close()

################################################################################################################################################
class frontend:
    # frontend class logger
    logging = create_logger("frontend")
    
    
    def __init__(self, feNoise = 51.735, adcResolution = 13, vref = 901.6e-3, lGain = 50.366, hGain = 503.66, chgShrCrossTalkMap = [0]):
        """
        Parameters
        fenoise         : {self.feNoise} 517.35 e
        
        """
        self.feNoise = feNoise
        self.adcResolution = adcResolution
        self.vref = vref
        self.lGain = lGain
        self.hGain = hGain
        self.chgShrCrossTalkMap = chgShrCrossTalkMap
        
        self.maxADC = 2**self.adcResolution - 1   
        self.adcScale = self.vref/self.maxADC #in units of V/count
        
        
        if 2*sum(chgShrCrossTalkMap) > 1.0:
            print("WARNING: Sum of charge sharing fractions cannot exceed 0.5. Cross talk will be disabled")
            self.chgShrCrossTalkMap = [0]    
            
            
        msg = f"""FRONT-END settings are:
        fenoise             : {self.feNoise} e
        adcresolution       : {self.adcResolution}-bit
        vref                : {self.vref} V
        lgain               : {self.lGain} mV/pC
        hgain               : {self.hGain} mV/pC
        chgShrCrossTalkMap  : {self.chgShrCrossTalkMap} %
        
        Calculated quantities:
        maxADC          : {self.maxADC}
        adcScale        : {self.adcScale}\n----------------
        """
        (self.logging).info(msg)
    
    # Apply the front-end noise to the profile of charge collected at the strips
    def applyFENoise(self, chgStripProfile: ROOT.TH1D, plotting=False) -> ROOT.TH1D:
        """
        Apply the front-end noise to the profile of charge collected at the strips
        
        Paramters
        ----------
            chgStripProfile (ROOT.TH1D) : profile of charge collected at the strips
            plotting (bool) : if True, plots the profile with the noise applied
            
        Returns
        -------
            chgStripPostProcessedProfile (ROOT.TH1D) : profile of charge collected at the strips with front-end noise applied
        """
        
        chgStripProfileCp = chgStripProfile.Clone()
        
        # Set the name, title and vertical axis title
        chgStripProfileCp.SetTitle(chgStripProfileCp.GetTitle().replace("Charge collected", "Charge collected (plus FE-noise)"))
        chgStripProfileCp.SetName(chgStripProfileCp.GetName().replace('Proj', 'ProjFE'))    #chgProjFEProfileY
        chgStripProfileCp.GetYaxis().SetTitle("charge collected with FE-noise [C]")
        chgStripProfileCp.GetXaxis().CenterTitle()
        chgStripProfileCp.GetYaxis().CenterTitle()
        
        # This array contains the front-end values only, with the purpose of be used for debugging/inspection
        if plotting: noiseVals = np.zeros(chgStripProfileCp.GetNbinsX())
        #print("Number of entries:", chgStripProfiCp.GetEntries())
        # Each bin content value is summed to a random value extracted from
        # a gaussian with mean 0 and sigma given by the feNoise parameter, converted in electric charge 
        for i in range(1, chgStripProfileCp.GetNbinsX() + 1):
            #print(f"Bin {i} has content {chgStripProfile.GetBinContent(i)}") 
            val = chgStripProfileCp.GetBinContent(i)
            newVal = val +  np.random.normal(0, self.feNoise * 1.602176e-19 )
            chgStripProfileCp.SetBinContent(i, newVal)
            
            if plotting: noiseVals[i-1] = newVal

            # Get existing error
            existing_error = chgStripProfileCp.GetBinError(i)
            
            # Calculate combined error
            combined_error = np.sqrt( existing_error**2 + (self.feNoise*1.602176e-19)**2 )
            # set bin error
            chgStripProfileCp.SetBinError(i, combined_error)

            # Debug 
            #print(newVal) 
            #print(val, np.random.normal(0, self.feNoise * 1.602176e-19 ))

        # Diagnostic plots
        if plotting == 10:
            view, ax = plt.subplots()
            view.suptitle("Profile with noise applied")
            ax.set_xlabel("strip no.")
            ax.set_ylabel("signal+noise [C]")
            ax.plot(noiseVals, label=f"noise set value is {np.round(self.feNoise,1)} e.")
            view.savefig("applyFENoise.pdf")
            np.save("thisPlot.bin", view, allow_pickle=True)
            plt.show()
            
        return chgStripProfileCp
    
    # TODO double check if the vaules sotred in .txt are correct
    def applyFENoisePerBunch(self, chgStripProfile: ROOT.TH1D, noise_range: np.array, plotting: bool) -> dict:
        """
        Apply the front-end noise to the profile of charge collected at the strips for each bunch
        
        Parameters
        ----------
            chgStripProfiles (dict) : dictionary of profiles of charge collected at the strips with bunch numbers as keys
            noise_range (tuple) : range of noise values to be applied, in the format (start, stop, step)
        
        Returns
        -------
            chgStripPostProcessedProfiles (dict) : dictionary of profiles of charge collected at the strips with front-end noise applied
        """

        average_noises = {}
        # Loop over the bunchID and noise values
        for bunchID, feNoise in enumerate(noise_range):
            # Initialize chgStripProfileCp using the input chgStripProfile
            chgStripProfileCp = chgStripProfile.Clone()

            # Set the name, title, and vertical axis title
            chgStripProfileCp.SetTitle(chgStripProfileCp.GetTitle().replace("Charge collected", f"Charge collected (plus FE-noise for Bunch {bunchID})"))
            chgStripProfileCp.SetName(chgStripProfileCp.GetName().replace('Proj', f'ProjFE_Bunch{bunchID}'))
            chgStripProfileCp.GetYaxis().SetTitle(f"charge collected with FE-noise for Bunch {bunchID} [C]")
            chgStripProfileCp.GetXaxis().CenterTitle()
            chgStripProfileCp.GetYaxis().CenterTitle()

            # This array contains the front-end values only, with the purpose of being used for debugging/inspection
            if (self.logging).level == 10:
                noiseVals = np.zeros(chgStripProfileCp.GetNbinsX())

            # Each bin content value is summed to a random value extracted from
            # a Gaussian with mean 0 and sigma given by the feNoise parameter, converted into electric charge
            for j in range(1, chgStripProfileCp.GetNbinsX() + 1):
                val = chgStripProfileCp.GetBinContent(j)
                rndmNoise = np.random.normal(0, feNoise * 1.602176e-19)
                newVal = val + rndmNoise
                chgStripProfileCp.SetBinContent(j, newVal)

                if (self.logging).level == 10:
                    noiseVals[j - 1] = newVal

                # Get existing error
                existing_error = chgStripProfileCp.GetBinError(j)

                # Calculate combined error
                combined_error = np.sqrt(existing_error ** 2 + (feNoise * 1.602176e-19) ** 2)
                # set bin error
                chgStripProfileCp.SetBinError(j, combined_error)
            
            # Diagnostic plots
            if (self.logging).level == 10:
                view, ax = plt.subplots()
                view2, bx = plt.subplots()
                view.suptitle("Profile with noise applied")
                view2.suptitle(f"Profile of noise applied for Bunch {bunchID}")
                ax.set_xlabel("strip no.")
                ax.set_ylabel("signal+noise [C]")
                bx.set_xlabel("strip no.")
                bx.set_ylabel("noise [C]")
                ax.plot(noiseVals, label=f"noise set value for Bunch {bunchID} is {np.round(feNoise, 1)} e.")
                bx.plot(rndmNoise, label=f"noise set value for Bunch {bunchID} is {np.round(feNoise, 1)} e.")
                view.show()
                view2.show()
                plt.show()
    
            average_noises[bunchID] = (rndmNoise, combined_error)

        # Write average noise values and errors to a .txt file
        with open('average_noises.txt', 'w') as f:
            for bunchID, (mean, std) in average_noises.items():
                f.write(f'Bunch {bunchID}:  noise = {mean}, error = {std}\n')
            f.close()

        return chgStripProfileCp


    # Apply the front-end amplification to the profile of charge projected at the strips
    def applyAmplification(self, chgStripWithFE: ROOT.TH1D) -> ROOT.TH1D:
        """
        Apply the front-end amplification to the profile of charge projected at the strips
        
        Paramters
        ----------
            chgStripWithFE (ROOT.TH1D) : profile of the charge projected at the strips with FE noise applied
            
        Returns
        -------
            lowGainVoltageProfile (ROOT.TH1D) : profile of low gain voltage projected at the strips with front-end low gain applied
            highGainVoltageProfile (ROOT.TH1D) : profile of high gain voltage projected at the strips with front-end high gain applied

        """
        lgADCinVoltProf, hgADCinVoltProf = chgStripWithFE.Clone(), chgStripWithFE.Clone()

        # Get sensor name, if upstream or downstream
        sensor = "upstream"
        if 'downstream' in chgStripWithFE.GetTitle(): sensor = "downstream"
        # Set title
        lgADCinVoltProf.SetTitle(f"Low gain ADC-input voltage {sensor} sensor")
        hgADCinVoltProf.SetTitle(f"High gain ADC-input voltage {sensor} sensor")
        # Set name (chgStripWithFE name is chgProjFEProfileY) is mapped into lgADCinVoltProfX or lgADCinVoltProfY
        lgADCinVoltProf.SetName("lgADCinVoltProf"+(chgStripWithFE.GetName())[-1])
        hgADCinVoltProf.SetName("hgADCinVoltProf"+(chgStripWithFE.GetName())[-1])
        # Set vertical axis title
        lgADCinVoltProf.GetYaxis().SetTitle("LG ADC-input (V)")
        hgADCinVoltProf.GetYaxis().SetTitle("HG ADC-input (V)")
        lgADCinVoltProf.GetXaxis().CenterTitle()
        lgADCinVoltProf.GetYaxis().CenterTitle()
        hgADCinVoltProf.GetXaxis().CenterTitle()
        hgADCinVoltProf.GetYaxis().CenterTitle()
        
        # Convert the charge to a voltage (the histograms contains chg. in C, the pre-factor 1e9 converts the mV/pC to V/C)
        lgADCinVoltProf.Scale(1e9*self.lGain)
        hgADCinVoltProf.Scale(1e9*self.hGain)

        # Diagnostic plots
        if (logging.level == 10):
            view, ax = plt.subplots()
            view.suptitle("Profile with amplification applied")
            ax.set_xlabel("strip no.")
            ax.set_ylabel("signal+noise [V]")
            ax.plot(lgADCinVoltProf, label=f"low gain set value is {np.round(self.lGain,1)} mV/pC.")
            ax.plot(hgADCinVoltProf, label=f"high gain set value is {np.round(self.hGain,1)} mV/pC.")
            view.legend(loc="upper right")
            view.show()
            plt.show()


        return (lgADCinVoltProf, hgADCinVoltProf)
    
    # TODO double check if the vaules sotred in .txt are correct
    # Apply the front-end amplification to the profile of charge projected at the strips
    def applyAmplificationPerBunch(self, chgStripWithFE: ROOT.TH1D, vRefRange: np.array) -> dict:
        """
        Apply the front-end amplification to the profile of charge projected at the strips for each bunch         
        Paramters
        ----------
            chgStripWithFE (ROOT.TH1D) : profile of the charge projected at the strips with FE noise applied
            vRefRange (np.array) : range of reference voltages to be applied
            
        Returns
        -------
            lowGainVoltageProfile (dict) : profile of low gain voltage projected at the strips with front-end low gain applied
            highGainVoltageProfile (dict) : profile of high gain voltage projected at the strips with front-end high gain applied

        """
        lgADCinVoltProf, hgADCinVoltProf = chgStripWithFE.Clone(), chgStripWithFE.Clone()

        # Get sensor name, if upstream or downstream
        sensor = "upstream"
        if 'downstream' in chgStripWithFE.GetTitle(): sensor = "downstream"
        # Set title
        lgADCinVoltProf.SetTitle(f"Low gain ADC-input voltage {sensor} sensor")
        hgADCinVoltProf.SetTitle(f"High gain ADC-input voltage {sensor} sensor")
        # Set name (chgStripWithFE name is chgProjFEProfileY) is mapped into lgADCinVoltProfX or lgADCinVoltProfY
        lgADCinVoltProf.SetName("lgADCinVoltProf"+(chgStripWithFE.GetName())[-1])
        hgADCinVoltProf.SetName("hgADCinVoltProf"+(chgStripWithFE.GetName())[-1])
        # Set vertical axis title
        lgADCinVoltProf.GetYaxis().SetTitle("LG ADC-input (V)")
        hgADCinVoltProf.GetYaxis().SetTitle("HG ADC-input (V)")
        lgADCinVoltProf.GetXaxis().CenterTitle()
        lgADCinVoltProf.GetYaxis().CenterTitle()
        hgADCinVoltProf.GetXaxis().CenterTitle()
        hgADCinVoltProf.GetYaxis().CenterTitle()
        
        # Convert the charge to a voltage (the histograms contains chg. in C, the pre-factor 1e9 converts the mV/pC to V/C)
        lgADCinVoltProf.Scale(1e9*self.lGain)
        hgADCinVoltProf.Scale(1e9*self.hGain)

        avglSigPerVref = {}
        avghSigPerVref = {}

        # Loop over the bunchID and noise values
        for bunchID, self.vref in enumerate(vRefRange):
            for i in range(1, lgADCinVoltProf.GetNbinsX() + 1):
                avglSigPerVref[bunchID] = (np.mean(lgADCinVoltProf.GetBinContent(i))/self.vref, np.std(lgADCinVoltProf.GetBinContent(i))/self.vref)
            for i in range(1, hgADCinVoltProf.GetNbinsX() + 1):
                avghSigPerVref[bunchID] = (np.mean(hgADCinVoltProf.GetBinContent(i))/self.vref, np.std(hgADCinVoltProf.GetBinContent(i))/self.vref)
                with open('average_signal_ADCRange.txt', 'w') as f:
                    for bunchID, (mean, std) in avglSigPerVref.items():
                        f.write(f'Bunch {bunchID}:  low gain signal = {mean}, error = {std}\n')
                    for bunchID, (mean, std) in avghSigPerVref.items():
                        f.write(f'Bunch {bunchID}:  high gain signal = {mean}, error = {std}\n')
                    f.close()

        # Diagnostic plots
        if (logging.level == 10):
            view, ax = plt.subplots()
            view.suptitle("Profile with amplification applied")
            ax.set_xlabel("strip no.")
            ax.set_ylabel("signal+noise [V]")
            ax.plot(lgADCinVoltProf, label=f"low gain set value is {np.round(self.lGain,1)} mV/pC.")
            ax.plot(hgADCinVoltProf, label=f"high gain set value is {np.round(self.hGain,1)} mV/pC.")
            view.legend(loc="upper right")
            view.show()
            plt.show()


        return (lgADCinVoltProf, hgADCinVoltProf)
    
    # Convert a profile from continuous voltage value to discrete value
    def applyADC(self, profile: ROOT.TH1D) -> ROOT.TH1D:
        """
        Apply the Analog-to-digital conversion to the profiles of low and high gain voltages at the strips.
        
        Paramters
        ----------
            contiProfile (ROOT.TH1D) : input 'continuous' profile (vertical unit is supposed to be V)
            
        Returns
        -------
            discrProfile (ROOT.TH1D) : output 'discrete' profile (vertical unit is ADC values)
        """
        
        # Debugging info about ADC module
        msg = f"""Applying AD conversion with parameters:
        adcResolution      : {self.adcResolution}-bit
        maxADC             : {self.maxADC} 
        vref               : {self.vref} V
        adcScale           : {self.adcScale} V/count\n----------------"""
        (self.logging).debug(msg)


        # Copy the input histogram
        out = profile.Clone()

        # Set the proper names, styles etc.
        out.SetTitle(out.GetTitle().replace("ADC-input voltage", "ADC-counts"))
        out.GetYaxis().SetTitle(f"ADC counts [0-{self.maxADC}]")
        out.SetName(out.GetName().replace('inVolt', 'cts'))
        
        # For each strip get the voltage after the amplification and convert it to ADC counts 
        # using the adcScale, finally, the value is clipped between 0 and the maximal allowed ADC value
        for i in range(out.GetNbinsX()):
            val = out.GetBinContent(i)
            newVal = np.clip(np.round(val / self.adcScale), 0, self.maxADC)
            out.SetBinContent(i, newVal)
        
        # Diagnostic plots
        if (logging.level == 10):
            view, ax = plt.subplots()
            view.suptitle("Profile of ADC counts")
            ax.set_xlabel("strip no.")
            ax.set_ylabel("ADC counts")
            ax.plot(out, label=f"ADC resolution is {self.adcResolution}-bit.")
            view.legend(loc="upper right")
            view.show()
            plt.show()
            
        return out
    
    # Attach the computed errors to the histograms
    def applyMCStatErrs(self, hist: dict, errs: dict)-> ROOT.TH1D:
        """
        Attach the computed errors to the histograms. The errors are computed from the MC statistics. 
        The function takes as input a histogram and an array of errors and sets the errors of the histogram to the values in the array.

        Parameters
        ----------
            hist (ROOT.TH1D) : input histogram
            errs (np.array) : array of errors to be attached to the histogram

        Returns
        -------
            hist (dict) : histogram with the errors attached
        """
        for bin in range(hist.GetNbinsX()):
            error = errs[bin]  # errs is a 1D array, so errs[bin-1] is a scalar
            if isinstance(error, np.ndarray) and error.size == 1:
                error = error.item()  # Convert size-1 numpy array to scalar
                hist.SetBinError(bin, error)
        return hist
          
          
          
    def digitizeRun(self, fname = "gammaRun_1M.root", bunchPartNb=100000)-> tuple:
        """
        Description
        -----------
        This function performs the digitization of the energy deposition maps. It takes as input the ROOT file with the energy deposition maps and returns a dictionary with the digitized histograms.

        Pipeline
        1. readEdepFromROOT - read the energy deposition from the ROOT file and store it in a dictionary
        2. computeError_v2 - compute the error on the charge profiles
        3. applyMCStatErrs - apply the errors to the histograms
        4. getChgProfiles - get the charge profiles from the energy deposition maps 
        5. applyFENoise - apply the front-end noise to the profile of charge collected at the strips
        6. applyAmplification - apply the front-end amplification to the profile of charge projected at the strips to get the input voltage
        7. applyADC - convert a voltage profile to a discrete ADC counts profile
        8. mapStripProfile - change labels to have strip numbers in the horizontal axis (instead of X/Y in mm)
        9. mapAllHistograms - map all histograms in the given dictionary using mapStripProfile
        10. dumpAllInRoot - dump all the objects in the args to the ROOT file with filename 'fname'
        11. simulateCrossTalk - simulate the cross talk effect between strips
        12. createTH1F - create a histogram for the stats
        13. createTGraph - create a TGraph

        Parameters
        ----------
            fname (str) : name of the ROOT file with the energy deposition maps
            bunchPartNb (int) : number of particles in the bunch

        Returns
        -------
            tuple : tuple of the digitized histograms
        """


        ## TODO Incorporate the computeError_v2 into the readEdepFromROOT(str, int) and applyMCStatErrs in getChgProfiles
        
        
        # Read the energy deposition from the ROOT file and store it in a dictionary
        bunchesData = readEdepFromROOT(fname, bunchPartNb)
        
        # Create dictionaries to store data for each bunch
        bunchIDs = []
        lgADCctsMeanUp = []
        lgADCctsStdDevUp = []
        hgADCctsMeanUp = [] 
        hgADCctsStdDevUp = []
        lgADCctsMeanDo = []
        lgADCctsStdDevDo = []
        hgADCctsMeanDo = [] 
        hgADCctsStdDevDo = []

        # Compute the error on the charge profiles
        bunchesData_err = computeError_v2(bunchesData)
        
        for bunchID, bunch_data in bunchesData.items():
            # Extract data from dictionary 
            edepMapUp = bunch_data[f"b{bunchID}_edepMapUp"]
            edepMapDo = bunch_data[f"b{bunchID}_edepMapDo"]

            # Get the profiles of charge deposited and projected
            chgProfilesUp, chgProfilesDo = getChgProfiles(edepMapUp), getChgProfiles(edepMapDo)
            
            # Attach the errors to the histograms
            for dir in [0,1]:
                self.applyMCStatErrs(chgProfilesUp['chgProjProfile'][dir], bunchesData_err[0, dir, :])
                self.applyMCStatErrs(chgProfilesDo['chgProjProfile'][dir], bunchesData_err[1, dir, :])


            # Apply strip cross-talk to the charge profiles
            #global _debug
            #_debug = True
            chgProfilesWithCT = simulateCrosstalk(chgProfilesUp, chgProfilesDo, self.chgShrCrossTalkMap)



            noise_range = np.linspace(51, 92, 10)
            # Apply the front-end noise to the profile of charge collected at the strips
            fEHistUp, fEHistDo = self.applyFENoise(chgProfilesUp['chgProjProfile'][0], noise_range, True), self.applyFENoise(chgProfilesDo['chgProjProfile'][1], noise_range, True)

            vRefRange = np.linspace(0.057, 0.504, 10)
            # Apply the front-end amplification to the profile of charge projected at the strips,
            lgHistUp, hgHistUp = self.applyAmplificationPerBunch(fEHistUp, vRefRange);  lgHistDo, hgHistDo = self.applyAmplificationPerBunch(fEHistDo, vRefRange)
            # Convert a profile from continuous voltage value to discrete value
            lgADCHistUp, hgADCHistUp = self.applyADC(lgHistUp), self.applyADC(hgHistUp); lgADCHistDo, hgADCHistDo = self.applyADC(lgHistDo), self.applyADC(hgHistDo)

            
            canvas1 = ROOT.TCanvas("canvas1", "getChgProfiles deposited")
            chgProfilesUp['chgDepProfile'][0].Draw("e1x0")


            canvas2 = ROOT.TCanvas("canvas2", "getChgProfiles projected")
            chgProfilesUp['chgProjProfile'][0].Draw("e1x0")
            
            
            canvas3 = ROOT.TCanvas("canvas3", "simulateCrossTalk")
            chgProfilesUp['chgProjWithCTProfile'][0].Draw("e1x0")
            
            
            canvas4 = ROOT.TCanvas("canvas4", "applyFENoisePerBunch")
            fEHistUp.Draw("e1x0")
            
            canvas5 = ROOT.TCanvas("canvas5", "applyAmplification")
            canvas5.Divide(2)
            lgHistUp.Draw("e1x0")
            canvas5.cd(2)
            lgADCHistUp.Draw("e1x0")
            
            input()
            exit()
            #observables = self.computeOnDigizedProfiles(lgADCHistUp, hgADCHistUp)
            
            #roFileClass.OPT_fill(bunch = 0, evt = 2, det=0, dir=0, optPar=np.array([0,1,2]), optParErr = np.array([0,1,2]))
            # Call the class to store these observables on the root file
            

            # Calculate stats
            bunchIDs.append(bunchID)
            # For upstream histograms
            lgADCctsMeanUp.append(mapStripProfile(lgADCHistUp).GetMean())
            lgADCctsStdDevUp.append(mapStripProfile(lgADCHistUp).GetStdDev())
            hgADCctsMeanUp.append(mapStripProfile(hgADCHistUp).GetMean())
            hgADCctsStdDevUp.append(mapStripProfile(hgADCHistUp).GetStdDev())
            # For downstream histograms
            lgADCctsMeanDo.append(mapStripProfile(lgADCHistDo).GetMean())
            lgADCctsStdDevDo.append(mapStripProfile(lgADCHistDo).GetStdDev())
            hgADCctsMeanDo.append(mapStripProfile(hgADCHistDo).GetMean())
            hgADCctsStdDevDo.append(mapStripProfile(hgADCHistDo).GetStdDev())

        return (edepMapUp, edepMapDo, chgProfilesUp, chgProfilesDo, lgHistUp, lgHistDo, hgHistUp, hgHistDo, fEHistUp, fEHistDo, lgADCHistUp, lgADCHistDo, hgADCHistUp, hgADCHistDo, bunchIDs, lgADCctsMeanUp, lgADCctsStdDevUp, hgADCctsMeanUp, hgADCctsStdDevUp, lgADCctsMeanDo, lgADCctsStdDevDo, hgADCctsMeanDo, hgADCctsStdDevDo)


    def computeOnDigizedProfiles(a, b):
        pass
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

if __name__=="__main__":
    #logging.setLevel(10)       # This correspond to debug mode
    
    
    # Define the front-end settings
    aFrontEnd = frontend(
        #feNoise = 51.735,
        adcResolution = 13,
        #vref = 901.6e-3,
        lGain = 50.366,
        hGain = 503.66,
        #Add crosstalk string here
        )
    
    # This is a container for the output of the digitization pipeline.
    # Meaning that you calculate the parameters you want to extract from the digitized profiles and you store them into this root file in such a way (RATIONAL AND LOGICAL WAY) that the
    # analysis that you have to perform in the (near) future is conveniently done using this root file.
    #roFileClass = rdataStruct_OPT("myDummyDigitizationOptimization.root")
    
    result = aFrontEnd.digitizeRun(bunchPartNb=10000, fname="build/dummyRun_100k.root")
