#################################################################################################
# @info Simple digitization for the LUXE GP                                                     #
# @date 23/10/26                                                                                #
#                                                                                               #
# This script provides a simple way to map the energy deposited in the GBP sensor to a          #
# statistics of digitized profiles where the charge transport, propagamtion, digitization and   #
# front-end effects are accounted for.                                                          #
# TRACE and VERBOSE levels used for very detailed debugging                                     #
#################################################################################################

from rdataStruct import rdataStruct_OPT
from multipledispatch import dispatch
from matplotlib import pyplot as plt
from tqdm import tqdm
import numpy as np
import ROOT
import os

from readFromMC import readFromMc
from frontend import frontend

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





# Change labels to have strip numbers in the horizontal axis (instead of X/Y in mm)
def mapStripProfile(hist: ROOT.TH1D) -> ROOT.TH1D:
    """
    TODO
    """




if __name__=="__main__":
    #logging.setLevel(10)       # This correspond to debug mode
    
    ## Read data from the MC
    #rMCClass = readFromMc("build/dummyRun_100k.root", 10000, 0.2, 27.0)
    #rMCClass.logging.setLevel(10)
    #rMCClass.readProjProfilesFromROOT()
    #exit()
    
    # Define the front-end settings
    aFrontEnd = frontend()
    result = aFrontEnd.digitizeRun(bunchPartNb=10000, fname="build/dummyRun_100k.root")
    
    

    
    # This is a container for the output of the digitization pipeline.
    # Meaning that you calculate the parameters you want to extract from the digitized profiles and you store them into this root file in such a way (RATIONAL AND LOGICAL WAY) that the
    # analysis that you have to perform in the (near) future is conveniently done using this root file.
    #roFileClass = rdataStruct_OPT("myDummyDigitizationOptimization.root")