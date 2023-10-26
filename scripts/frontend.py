#################################################################################################
# @info description                                                     #
# @date 23/10/26                                                                                #
#                                                                                               #

#################################################################################################
from multipledispatch import dispatch
from matplotlib import pyplot as plt
from logger import create_logger
import numpy as np
import ROOT




class frontend:
    # frontend class logger
    logging = create_logger("frontend")
    
    
    @dispatch(float, int, float, float, float, list)
    def __init__(self, feNoise: float, adcResolution: int, vref: float, lGain: float, hGain: float, chgShrCrossTalkMap: list):
        """
        Create a front-end class representing the effect of a certain front-end setup (electronics) to the initial profile
        
        Parameters
        ----------
            feNoise (float) : Front-end electronics noise (in electron units)
            adcResolution (float) : Bit depth of the ADC converter (number of bits)
            vref (float) : Delete me
            lGain (float) : Pre-amplifier gain for low gain channel (in mV/pC)
            hGain (float) : Pre-amplifier gain for high gain channel (in mV/pC)
            chgShrCrossTalkMap (float) : List with percentages of the cross-talk charge sharing between the nearest neighbours strips (es. [0.1,0.002, 0.0003] means that 0.1 is shared between strips at 1 distance, 0.002 between strips at distance 2 and so on)
        
        """
        
        # Set internal variables of the class with the external parameters
        self.feNoise = feNoise
        self.adcResolution = adcResolution
        self.vref = vref
        self.lGain = lGain
        self.hGain = hGain
        self.chgShrCrossTalkMap = chgShrCrossTalkMap
        
        self.maxADC = 2**self.adcResolution - 1   
        self.adcScale = self.vref/self.maxADC #in units of V/count
        
        
        if 2*sum(chgShrCrossTalkMap) > 1.0:
            (self.logging).warning("Sum of charge sharing fractions cannot exceed 0.5. Cross talk will be disabled")
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
    
    # Overload of the init class function for taking an input projProfile profile and processing with default parameters
    @dispatch(ROOT.TH1D, float, int, float, float, float, list)
    def __init__(self, projChgProfile: ROOT.TH1D, feNoise: float, adcResolution: int, vref: float, lGain: float, hGain: float, chgShrCrossTalkMap: list):
        (self.logging).error("I'm not implemented yet")
        raise Exception("I'm not implemented yet")
    
    
    # Overload of the init class function for loading profile from scratch and with default parameters
    @dispatch()
    def __init__(self):
        return (self.__init__)(feNoise = 51.735, adcResolution = 13, vref = 901.6e-3, lGain = 50.366, hGain = 503.66, chgShrCrossTalkMap = [0])
    
    
    
    
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
        if ((self.logging).level == 10):
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
        if ((self.logging).level == 10):
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
        if ((self.logging).level == 10):
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
    
    
    
    # Simulates the cross talk between two strips of the sensors. The cross talk is simulated by adding a normal distribution centered around the selected strip number for each sensor.
    @dispatch(ROOT.TH1D, list)
    def simulateCrosstalk(self, chgProfile: ROOT.TH1D, chgShrNN: list) -> ROOT.TH1D:
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
        if ((self.logging).level == 10):
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
        return chgProfileWithCrossTalk

    # Overload to accomodate the structure of the other functions
    @dispatch(dict, dict, list)
    def simulateCrosstalk(self, chgProfilesUp: dict, chgProfilesDo: dict, chgShrNN: list) -> tuple:
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
        chgProjWithCTProfileUp_X = self.simulateCrosstalk(chgProfilesUp['chgProjProfile'][0], chgShrNN)
        chgProjWithCTProfileUp_Y = self.simulateCrosstalk(chgProfilesUp['chgProjProfile'][1], chgShrNN)
        chgProfilesUp['chgProjWithCTProfile'] = (chgProjWithCTProfileUp_X, chgProjWithCTProfileUp_Y)
        
        chgProjWithCTProfileDo_X = self.simulateCrosstalk(chgProfilesUp['chgProjProfile'][0], chgShrNN)
        chgProjWithCTProfileDo_Y = self.simulateCrosstalk(chgProfilesUp['chgProjProfile'][1], chgShrNN)
        chgProfilesDo['chgProjWithCTProfile'] = (chgProjWithCTProfileDo_X, chgProjWithCTProfileDo_Y)
        
        return (chgProfilesUp, chgProfilesDo)
    
    
    
          
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
        bunchesData_err = calculateProfileStatErrs(bunchesData)
        
        for bunchID, bunch_data in bunchesData.items():
            # Extract data from dictionary 
            edepMapUp = bunch_data[f"b{bunchID}_edepMapUp"]
            edepMapDo = bunch_data[f"b{bunchID}_edepMapDo"]

            # Get the profiles of charge deposited and projected
            chgProfilesUp, chgProfilesDo = getChgDepProjProfiles(edepMapUp), getChgDepProjProfiles(edepMapDo)
        
            # Attach the errors to the histograms
            for dir in [0,1]:
                self.applyMCStatErrs(chgProfilesUp['chgProjProfile'][dir], bunchesData_err[0, dir, :])
                self.applyMCStatErrs(chgProfilesDo['chgProjProfile'][dir], bunchesData_err[1, dir, :])


            # Apply strip cross-talk to the charge profiles
            #global _debug
            #_debug = True
            chgProfilesWithCT = self.simulateCrosstalk(chgProfilesUp, chgProfilesDo, self.chgShrCrossTalkMap)



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