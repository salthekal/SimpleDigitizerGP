#################################################################################################
# @info Simple digitization for the LUXE GP                                                     #
# @date 23/10/26                                                                                #
#                                                                                               #
# This script provides a simple way to map the energy deposited in the GBP sensor to a          #
# statistics of digitized profiles where the charge transport, propagamtion, digitization and   #
# front-end effects are accounted for.                                                          #
# TRACE and VERBOSE levels used for very detailed debugging                                     #
#################################################################################################
from multipledispatch import dispatch
from matplotlib import pyplot as plt
from tqdm import tqdm
import numpy as np
import ROOT
import os

from readFromMC import readFromMc
from frontend import frontend
from featureExtractor import featureExtractor

# digiManager.py logger
from logger import create_logger
logging = create_logger("digiManager")
    
# Enable ROOT implicit multithreading
ROOT.EnableImplicitMT(16)


"""
TODO
1. Integrate the rdataStruct_OPT class
2. Memory management for ROOT objects to cleanup useless stuffs
3. Multithreading for parallel processing of the parameter space for the readMC and frontend classes
4. Setup the jobs for the calculation
5. Write down a class which calculated on the digitized profiles the parameters that we want to extract from the digitization pipeline
"""

# Converts the X/Y position to a strip number between 1-192. By convention, strip #1 is the one with negative X/Y coordinate
def calculateStripNb(pos_mm: float, stripPitch_mm: float, nbStrips: int) -> int:
    """
    Converts the X/Y position to a strip number between 1-192. By convention, strip #1 is the one with negative X/Y coordinate
    
    Parameters
    ----------
        pos_mm (float) : position in mm
        stripPitch_mm (float) : strip pitch in mm
        nbStrips (int) : total number of strips in the sensor
        
    Returns
    -------
        stripNb (int) : strip number between 1-192
    """
    
    return int( (pos_mm + stripPitch_mm*nbStrips/2.0) / (stripPitch_mm*nbStrips) * nbStrips + 1.0 )


@dispatch(str, str, int, float, float, float, int, float, float, list)
def pipeline(mcPath: str, outPath: str, _bunchParNb: int, _cce: float, _avgPairEn: float, _fNoise: float, _iADCres: int, _fOlScale: float, _fGain: float, _vChgShrCrossTalkMap: list):
    """
    Apply the full pipeline from the incident radiation (readFromMc) to the digitized signal (ADC)
    
    Parameters
    ----------
        mcPath (str) : Filename of the MC ROOT file containing the data
        outPath (str) : Filename of the output ROOT file where the parameters from the feature extraction are stored
        _bunchParNb (int) : Number of particles in a bunch
        _cce (float) : Charge collection efficiency of the sensor for the bare geometrical projection of the dep chg. to proj chg.
        _avgPairEn (float) : Average energy required to create an electron-hole pair in the sensor (in eV)
        _fNoise (float) : Front-end electronics noise (in electron units) and internally converted in Coulomb
        _iADCres (int) : Number of bits of the analog-to-digital converter (12-bit -> 12)
        _fOlScale (float) : Full scale range that the frontend is capable of digitizing (in electron units)
        _fGain (float) : Amplifier gain (ration between final charge and initial charge)
        _vChgShrCrossTalkMap (list) : Cross-talk charge sharing percentages (es. [0.1,0.002, 0.0003] means that 0.1 is shared between strips at 1 distance, 0.002 between strips at distance 2, etc.)
        
    Returns
    -------
        None     
    """
    
    # Create the pipeline parameters
    _pipelinePars = {
        "bunchParNb" : _bunchParNb,
        "cce" : _cce,
        "avgPairEn" : _avgPairEn,
        "fNoise" : _fNoise,
        "iADCres" : _iADCres,
        "fOlScale" : _fOlScale,
        "fGain" : _fGain,
        "vChgShrCrossTalkMap" : _vChgShrCrossTalkMap
    }
    
    # Read data from the MC
    rMCClass = readFromMc(rFname = mcPath, bunchParNb = _bunchParNb, cce = _cce, avgPairEn = _avgPairEn)
    projChgProfs = rMCClass.readProjProfilesFromROOT()

    # Define the front-end settings
    #aFrontEnd = frontend(fNoise = (1e-15/1.60e-19), iADCres = 13, fOlScale = (12e-15/1.60e-19), fGain = 1.0, vChgShrCrossTalkMap = [0], projChgProfiles = tmp)
    aFrontEnd = frontend(projChgProfiles = projChgProfs, fNoise = _fNoise, iADCres = _iADCres, fOlScale = _fOlScale, fGain = _fGain, vChgShrCrossTalkMap = _vChgShrCrossTalkMap)
    aFrontEnd.doPipeline()
    
    # The following class calculated the observables that one is interested to extract from the profiles
    fExInstance = featureExtractor(roFname = outPath, dgtChgProfiles = aFrontEnd.projChgProfiles, pipelinePars = _pipelinePars)
    fExInstance.fitSchemeA()
    
    fExInstance.writeFeatures()
    logging.info("Done")


# Overload of the pipeline function with default parameters
@dispatch()
def pipeline():
    """
    Overload of the pipeline function with default parameters
    
    Parameters
    ----------
    1. mcPath = "build/dummyRun_100k.root",
    2. outPath = "testFeatureExtractor.root"
    3. _bunchParNb = 10000,
    4. _cce = 0.2,
    5. _avgPairEn = 27.0,
    6. _fNoise = 0,
    7. _iADCres = 8, 
    8. _fOlScale = (10e-15/1.60e-19),
    9. _fGain = 1.0,
    10. _vChgShrCrossTalkMap = [0]
    """
    pipeline(
        "build/dummyRun_100k.root",
        "testFeatureExtractor.root",
        10000,
        0.2,
        27.0,
        0.0,
        8, 
        (10e-15/1.60e-19),
        1.0,
        [0])







if __name__=="__main__":
    #logging.setLevel(10)       # This correspond to debug mode
    
    pipeline()
    exit()
     
    # Canvas to plot the digitized profile
    canvas1 = ROOT.TCanvas("canvas1", "canvas1")
    projChgProfs[0]['d0_x'].Draw("AP")
    canvas1.Update()
    
    canvas2 = ROOT.TCanvas("canvas2", "canvas2")
    aFrontEnd.projChgProfiles[0]['d0_x'].Draw("AP")
    canvas2.Update()
    input()
    exit()
        
