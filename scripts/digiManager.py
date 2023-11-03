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
from sys import argv as CLIargs
from inspect import signature
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
        "vChgShrCrossTalkMap" : np.array(_vChgShrCrossTalkMap + [0]*(10-len(_vChgShrCrossTalkMap)))
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
        [0.1,0.01,0.001])



# Generate a parameter space. The idea of this function is to optimize the phase space w.r.t. the naive linspace tensor product
def makePhaseSpace(cce: tuple, avgPairEn: tuple, fNoise: tuple, iADCres: tuple, fOlScale: tuple, fGain: tuple, vChgShrCrossTalkMap: tuple) -> np.array:
    """
    Expand the input tuples and generate the phase space to sample with the simulation
    
    Parameters
    ----------
        cce (tuple) : Range for the cce in the form (cceStart, cceEnd, Nppoints). (It includes the endpoint)
        avgPairEn (tuple) : Range for the avgPairEn in the form (avgPairEnStart, avgPairEnEnd, Nppoints). (It includes the endpoint)
        fNoise (tuple) : Range for the fNoise in the form (fNoiseStart, fNoiseEnd, Nppoints). (It includes the endpoint)
        iADCres (tuple) : Range for the iADCres in the form (iADCresStart, iADCresEnd, Nppoints). (It includes the endpoint)
        fOlScale (tuple) : Range for the fOlScale in the form (fOlScaleStart, fOlScaleEnd, Nppoints). (It includes the endpoint)
        fGain (tuple) : Range for the fGain in the form (fGainStart, fGainEnd, Nppoints). (It includes the endpoint)
        vChgShrCrossTalkMap (tuple) : Range for the vChgShrCrossTalkMap in the form (vChgShrCrossTalkMapStart, vChgShrCrossTalkMapEnd, Nppoints). (It includes the endpoint)
    
    Returns
    -------
        phaseSpace (np.array) : stack with the phase space linstack
    """
    
    # Expand the input tuples
    cceA, cceB, cceN = cce
    avgPairEnA, avgPairEnB, avgPairEnN = avgPairEn
    fNoiseA, fNoiseB, fNoiseN = fNoise
    iADCresA, iADCresB, iADCresN = iADCres
    fOlScaleA, fOlScaleB, fOlScaleN = fOlScale
    fGainA, fGainB, fGainN = fGain
    vChgShrCrossTalkMapA, vChgShrCrossTalkMapB, vChgShrCrossTalkMapN = vChgShrCrossTalkMap
    
    # Throw exception if the parameters are unallowed
    if len(vChgShrCrossTalkMapA)!=len(vChgShrCrossTalkMapB):
        logging.critical(f"The length of vChgShrCrossTalkMapA and vChgShrCrossTalkMapB should be equal. Instead you have used {vChgShrCrossTalkMapA} and {vChgShrCrossTalkMapB} with length {len(vChgShrCrossTalkMapA)} and {len(vChgShrCrossTalkMapB)}, respectively.")
        raise Exception("Invalid vChgShrCrossTalkMapA or vChgShrCrossTalkMapB (len).")
    
    # Throw some warnings if the picked phase space is idiotic
    
    # Generate the linear spaces in a sensible way
    phaseSpace_cce = np.linspace(cceA, cceB, cceN, endpoint=True)
    phaseSpace_avgPairEn = np.linspace(avgPairEnA, avgPairEnB, avgPairEnN, endpoint=True)
    phaseSpace_fNoise = np.linspace(fNoiseA, fNoiseB, fNoiseN, endpoint=True)
    phaseSpace_iADCres = np.linspace(iADCresA, iADCresB, iADCresN, endpoint=True)
    phaseSpace_fOlScale = np.linspace(fOlScaleA, fOlScaleB, fOlScaleN, endpoint=True)
    phaseSpace_fGain = np.linspace(fGainA, fGainB, fGainN, endpoint=True)
    phaseSpace_vChgShrCrossTalkMap = np.linspace(vChgShrCrossTalkMapA, vChgShrCrossTalkMapB, vChgShrCrossTalkMapN, endpoint=True)
    
    return np.stack(phaseSpace_cce, phaseSpace_avgPairEn, phaseSpace_fNoise, phaseSpace_iADCres, phaseSpace_fOlScale, phaseSpace_fGain, phaseSpace_vChgShrCrossTalkMap)
    
    


def makeJobs():
    pass



if __name__=="__main__":
    pass
    #logging.setLevel(10)       # This correspond to debug mode
    #pipeline()
    #exit()
    





########################################################################################################################
########################################################################################################################
# Documentation
def printHelp():
    """
    Print the CLI documentation
    """
    # Take a copy of the list of symbols defined in the global scope
    globalSymbols = dict(globals())
    
    print("Usage: python digiManager.py <function_name> <arg1> <arg2> ...")
    print("Available functions are: ")
    count = 1
    for symbolName in globalSymbols:
        if symbolName == "printHelp": continue
        # Check if the symbol is a callable
        if callable(globals()[symbolName]):
            if count > 4: print(f"\t{count-4}. {symbolName}")
            count += 1

##############################
# Handle CLI arguments
if len(CLIargs) < 2:
    # Documentation
    printHelp()
else:
    # Take the list of symbols defined in the global scope
    globalSymbols = globals()

    # Parse the first argument
    callableName = CLIargs[1]
    if (callableName not in globalSymbols):
        logging.error(f"Symbol {callableName} not found in the global scope")
        printHelp()
        exit(-1)
    elif not callable(globals()[callableName]):
        logging.error(f"Symbol {callableName} in globals but not a callable")
        printHelp()
        exit(-1)
    
    
    # Get the function object
    function_to_call = globals()[callableName]
    
    def get_signatures(dispatched_function):
        """
        Get the list of signatures for the dispatched function. Generalization of the signature function for the multidispatch
        """
        # Get the Signature of the callable and the input parameters
        function_signature = signature(dispatched_function)
        
        if (f"{function_signature}" == "(*args, **kwargs)"):
            # dispatch
            signatures = [entry[0] for entry in dispatched_function.funcs.items()]
        else:
            signatures = [entry[1].annotation for entry in function_signature.parameters.items()]
            signatures = list([tuple(signatures)])
        return signatures
    
    
    # Get the Signature of the callable and the input parameters and extract the names and types of parameters
    param_info = get_signatures(function_to_call)
    
    # Extract the function arguments from the command line
    args_from_cli = CLIargs[2:]
    
    # Check if the number of arguments matches the number of parameters
    argNbs_NotMatching = True
    for entry in param_info:
        if len(entry) == len(args_from_cli):
            argNbs_NotMatching = False
            break
    if argNbs_NotMatching:
        logging.error("Number of arguments does not match the function signature.")
        print(param_info)
        exit(-1)


    # Attempt to convert arguments to the expected types
    ## Generalization for multipledispatch
    argConvertedFail = True
    for par_types in param_info:
        toCastParNb = 0
        converted_args = []
        for argI, par_type in enumerate(par_types):
            #logging.debug(par_type, argI, args_from_cli[argI])
            try:
                if args_from_cli[argI].isdecimal() and par_type is str: continue
                if (args_from_cli[argI].replace('.', '').isdecimal()) and par_type in [float, np.double, np.single]: continue
                castedPar = par_type(args_from_cli[argI])
                converted_args.append(castedPar)
                toCastParNb += 1
            except:
                #logging.debug(f"Failed to convert {args_from_cli[argI]} to {par_type}. breaking to next set of args")
                break
        if len(converted_args) == len(args_from_cli):
            argConvertedFail = False
            break        
    if argConvertedFail:
        logging.error("Unable to convert the arguments to any of the expected types.")
        exit(-1)
        

    # Info message
    msg = f"Running: {callableName} with arguments "
    for arg in converted_args: msg += f"{arg} "
    msg = msg[:-1]
    logging.info(msg)
    
    # Call the function with the converted arguments
    function_to_call(*converted_args)
    
    # Goodbye
    msg = ["Whatever happens, happens. - Spike Spiegel", "Everything has a beginning and an end. - Jet Black" , "They say hunger is the best spice. - Spike Spiegel", "That's what she said. - M. Scott"]
    from random import choice
    logging.status(choice(msg)+". Goodbye :) \n")

########################################################################################################################
########################################################################################################################