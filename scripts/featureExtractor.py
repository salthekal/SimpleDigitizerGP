#################################################################################################
# @info Extract the relevant features from the digitized profiles of charge deposited in the    #
#       upstream and downstream sensors                                                         #
# @date 23/10/31                                                                                #
#                                                                                               #
# The list of features are for now limited to the very basic set of parameters that you read    #
# by the (gaussian + constant) fit parameters of the digitized charge profiles.                 #
#################################################################################################
from rdataStruct import rdataStruct_OPT
from multipledispatch import dispatch
from logger import create_logger
import numpy as np
import ROOT
import os

# This is a container for the output of the digitization pipeline.
# Meaning that you calculate the parameters you want to extract from the digitized profiles and you store them into this root file in such a way (RATIONAL AND LOGICAL WAY) that the
# analysis that you have to perform in the (near) future is conveniently done using this root file.
#roFileClass = rdataStruct_OPT("myDummyDigitizationOptimization.root")


class featureExtractor:
    """
    Extended class description
    """
    # frontend class logger
    logging = create_logger("frontend")
    
    
    def initialize(self, roFname: str):
        """
        Define the gaussian fit and other ROOT based objects that should be persistent in the class memory space
        """
        
        # Define the fit schemes
        ## Scheme A
        self.fitSchemeA_gaussFit = ROOT.TF1("fitSchemeA", "gaus(0) + pol0(3)", -10.0, 10.0)
        self.fitSchemeA_gaussFit.SetParLimits(0,   0., 8191.)        # Constant
        self.fitSchemeA_gaussFit.SetParLimits(1, -10.,   10.)        # Mean
        self.fitSchemeA_gaussFit.SetParLimits(2,0.100, 10./3)        # Sigma
        self.fitSchemeA_gaussFit.SetParLimits(3,   0.,  100.)        # Baseline
        self.fitSchemeA_gaussFit.SetParNames("fSA_amp", "fSA_mea", "fSA_sig", "fSA_bck")
        #self.fitSchemeA_gaussFit.SetParameters(8191./2., 92., 2.0, 0.)
        
        # Create the class containing the output of the feature extraction data
        self.roClass = rdataStruct_OPT(roFname)
           
    
    def __init__(self, dgtChgProfiles: np.array, roFname: str, pipelinePars: dict):
        """
        Description
        
        Parameters
        ----------
            dgtChgProfiles (np.array) : list of digitized charge profiles for the upstream and downstream detector in the horizontal and vertical direction, respectively
            roFname (str) : Filename of the output ROOT file where the parameters from the feature extraction are stored
            pipelinePars (dict) : Dictionary with the parameters used in the digitization pipeline 

        """
        
        # Copy the list of TGraphErrors internally so that the input objects are not modified
        self.dgtChgProfiles = np.array([{'d0_x' : item['d0_x'].Clone(), 'd1_y': item['d1_y'].Clone()} for item in dgtChgProfiles])
        
        # if the roFname already exists or the filename is invalid, then default to the dgtOut.root in the current dir
        if os.path.exists(roFname):
            (self.logging).warning(f"File {roFname} already exists. Defaulting to dgtOut.root to prevent overwrite.")
            roFname = "dgtOut.root"
        elif roFname.find('.root') == -1:
            (self.logging).error(f"Filename {roFname} is invalid. Defaulting to dgtOut.root to prevent overwrite.")
            roFname = "dgtOut.root"
        # Call the initialize function
        self.initialize(roFname)
        
        # If the class if called without any indication of the paraemters used in the digitization pipeline, throw a warning
        if len(pipelinePars) == 0: (self.logging).error("No information about digitization pipeline parameters. Check!")
        self.dgtPars = pipelinePars
        
        # Logging
        msg = f"""Class featureExtractor will look for:
        1. Profile characteristics from a (gaussian+pol0) fit (amplitude, mean, sigma, constant)
        2. ?\n----------------
        """
        (self.logging).info(msg)
    
    
    # Initialize and perform the gaussian fit
    @dispatch(ROOT.TGraphErrors)
    def fitSchemeA(self, profile: ROOT.TGraphErrors) -> dict:
        """
        Fit the profile using the fit scheme called self.fitSchemeA_gaussFit
        
        Parameters
        ----------
            profile (ROOT.TGraphErrors) : profile to fit
        
        Returns
        -------
            fitPars (dict) : dictionary with fit parameters and their errors
        """
        
        # Copy the fitSchemeA into a new object (otherwise the ParLimits and other attributes of the fit will be changed)
        fitSchemeA_gaussFit = self.fitSchemeA_gaussFit.Clone()
        
        # Calculate the optimal values to correctly inizialize the fit parameters
        ## Mean and standard deviation
        grPoints =np.column_stack((profile.GetX(), profile.GetY()))
        grMin, grMax = np.min(grPoints[:,1]), np.max(grPoints[:, 1])
        grStdX, grStdY = profile.GetRMS(1), profile.GetRMS(2)
        grMeanX, grMeanY = np.sum(grPoints[:, 0] * grPoints[:, 1]) / np.sum(grPoints[:, 1]), profile.GetMean(2)
        
        #print("grMin", grMin, grMax)
        #print("grStdX", grStdX, grStdY)
        #print("grMeanX", grMeanX, grMeanY)
        #print("-------")
        
        ## Calculate the intervals for the parameters
        amplL, amplH = grMin,               8192.0
        meanL, meanH = (grMeanX-grStdX),    (grMeanX+grStdX)
        sigmL, sigmH = 0.01*grStdX,         2.0*grStdX
        baseL, baseH = 0,                   1.2*grMin
        
        #print(amplL, amplH)
        #print(meanL, meanH)
        #print(sigmL, sigmH)
        #print(baseL, baseH)
        
        # Set the optimal values to correclty initialize the fit
        fitSchemeA_gaussFit.SetRange(grMeanX-2*grStdX, grMeanX+2*grStdX)
        fitSchemeA_gaussFit.SetParLimits(0, amplL, amplH)       # Constant
        fitSchemeA_gaussFit.SetParLimits(1, meanL, meanH)       # Mean
        fitSchemeA_gaussFit.SetParLimits(2, sigmL, sigmH)       # Sigma
        fitSchemeA_gaussFit.SetParLimits(3, baseL, baseH)       # Baseline
        fitSchemeA_gaussFit.SetParameters( 0.5*(amplL+amplH), 0.5*(meanL+meanH), 0.5*(sigmL+sigmH), 0.5*(baseL+baseH))
        
        # Fit the profile
        fitOpts  = "Q  "   # Quiet mode
        fitOpts += "R "   # Use the range specified in the function range
        fitOpts += "S "   # The result of the fit is returned in the TFitResultPtr
        #profile.Draw("AP")
        fitResultPtr = profile.Fit(fitSchemeA_gaussFit, fitOpts)
        # Get the covariance matrix, from which the chi2 and fit pars can be read 
        #covMatrix = fitResultPtr.GetCovarianceMatrix()
        
        # Chi2
        fitPars = {}
        fitPars['chi2'] = fitSchemeA_gaussFit.GetChisquare()
        # Degrees of freedom
        fitPars['ndf'] = fitSchemeA_gaussFit.GetNDF()
        # Reduced chisquare
        fitPars['rchi2'] = fitPars['chi2'] / fitPars['ndf']
        # Get fit parameters and errors in the fit result
        for i in range(4):
            parName = fitSchemeA_gaussFit.GetParName(i)
            fitPars[parName] = (fitSchemeA_gaussFit.GetParameter(i), fitSchemeA_gaussFit.GetParError(i))
        
        return fitPars
    

    # Overload of the fitSchemeA function
    @dispatch()
    def fitSchemeA(self):
        for entry in self.dgtChgProfiles:
            entry['d0_fitSchemeA'] = self.fitSchemeA(entry['d0_x'])
            entry['d1_fitSchemeA'] = self.fitSchemeA(entry['d1_y'])
            
    # Write the results of the feature extraction into the root file
    def writeFeatures(self):
        for i, entry in enumerate(self.dgtChgProfiles):
            data = {}
            data['bunch'] = i
            data['evt'] = i
            data['det'] = np.array([0,1], dtype=np.uint32)
            
            data.update(self.dgtPars)
            
            data['chi2'] = np.array([entry['d0_fitSchemeA']['chi2'], entry['d1_fitSchemeA']['chi2']])
            data['ndf'] = np.array([entry['d0_fitSchemeA']['ndf'], entry['d1_fitSchemeA']['ndf']])
            data['rchi2'] = np.array([entry['d0_fitSchemeA']['rchi2'], entry['d1_fitSchemeA']['rchi2']])
            data['fSA_amp'] = np.array([entry['d0_fitSchemeA']['fSA_amp'][0], entry['d1_fitSchemeA']['fSA_amp'][0]])
            data['fSA_mea'] = np.array([entry['d0_fitSchemeA']['fSA_mea'][0], entry['d1_fitSchemeA']['fSA_mea'][0]])
            data['fSA_sig'] = np.array([entry['d0_fitSchemeA']['fSA_sig'][0], entry['d1_fitSchemeA']['fSA_sig'][0]])
            data['fSA_bck'] = np.array([entry['d0_fitSchemeA']['fSA_bck'][0], entry['d1_fitSchemeA']['fSA_bck'][0]])
            data['fSA_amp_err'] = np.array([entry['d0_fitSchemeA']['fSA_amp'][1], entry['d1_fitSchemeA']['fSA_amp'][1]])
            data['fSA_mea_err'] = np.array([entry['d0_fitSchemeA']['fSA_mea'][1], entry['d1_fitSchemeA']['fSA_mea'][1]])
            data['fSA_sig_err'] = np.array([entry['d0_fitSchemeA']['fSA_sig'][1], entry['d1_fitSchemeA']['fSA_sig'][1]])
            data['fSA_bck_err'] = np.array([entry['d0_fitSchemeA']['fSA_bck'][1], entry['d1_fitSchemeA']['fSA_bck'][1]])
            
            (self.roClass).OPT_fill(**data)
        
        (self.roClass).closeROOT()