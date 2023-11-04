// Include the necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

// Function to create a histogram from the TTree
void CreateHistogram(TTree* tree, const char* branchName, const char* condition, const char* drawOption, const char* title, const char* xTitle, const char* yTitle, const char* saveAs) {
    // Create a canvas to draw the histogram
    TCanvas* canvas = new TCanvas();
    
    // Draw the histogram
    tree->Draw(Form("%s>>hist", branchName), condition, drawOption);
    
    // Get the histogram from the current directory
    TH1D* hist = (TH1D*)gDirectory->Get("hist");
    
    // Set the titles
    hist->SetTitle(title);
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);
    
    // Save the histogram
    canvas->SaveAs(saveAs);
    
    // Clean up
    delete canvas;
    delete hist;
}

void CreateHistograms() {
    // Open the .root file
    TFile* file = TFile::Open("dgtOut.root");
    
    // Get the TTree
    TTree* tree = (TTree*)file->Get("OPT");
    
    // Create histograms
    CreateHistogram(tree, "fNoise", "(det==0)", "hist", "Noise in the Upstream sensor Histogram", "Noise", "Counts", "fNoiseU.png");
    CreateHistogram(tree, "fNoise", "(det==1)", "hist", "Noise in the Downstream sensor Histogram", "Noise", "Counts", "fNoiseD.png");

    CreateHistogram(tree, "fNoise", "(det==0)*bunch", "hist", "Noise per Bunch in the Upstream Histogram", "Bunch ID", "Noise", "fNoiseBunchU.png");
    CreateHistogram(tree, "fNoise", "(det==1)*bunch", "hist", "Noise per Bunch in the Downstream Histogram", "Bunch ID", "Noise", "fNoiseBunchD.png");

    CreateHistogram(tree, "fOlScale", "(det==0)", "hist", "Fullscale Range in the Upstream Histogram", "Fullscale Range", "Counts", "fOlScaleU.png");
    CreateHistogram(tree, "fOlScale", "(det==1)", "hist", "Fullscale Range in the Downstream Histogram", "Fullscale Range", "Counts", "fOlScaleD.png");

    CreateHistogram(tree, "fOlScale", "(det==0)*bunch", "hist", "Fullscale Range per Bunch in the Upstream Histogram", "Fullscale Range Bunch", "Counts", "fOlScaleBunchU.png");
    CreateHistogram(tree, "fOlScale", "(det==1)*bunch", "hist", "Fullscale Range per Bunch in the Downstream Histogram", "Fullscale Range Bunch", "Counts", "fOlScaleBunchD.png");


    CreateHistogram(tree, "fGain", "(det==0)", "hist", "Gain in the Upstream Histogram", "Gain", "Counts", "fGainU.png");
    CreateHistogram(tree, "fGain", "(det==1)", "hist", "Gain in the Downstream Histogram", "Gain", "Counts", "fGainD.png");


    CreateHistogram(tree, "fGain", "(det==0)*bunch", "hist", "Gain Bunch in the Upstream Histogram", "Gain Bunch", "Counts", "fGainBunchU.png");
    CreateHistogram(tree, "fGain", "(det==1)*bunch", "hist", "Gain Bunch in the Downstream Histogram", "Gain Bunch", "Counts", "fGainBunchD.png");


    // Observables

    CreateHistogram(tree, "fNoise", "(det==0)*(fSA_sig_err)", "hist", "Noise level vs sigma error in the Upstream sensor ", "Noise level", "Sigma error", "fNoiseSigErrorU.png");
    CreateHistogram(tree, "fNoise", "(det==1)*(fSA_sig_err)", "hist", "Noise level vs sigma error in the Downstream sensor ", "Noise level", "Sigma error", "fNoiseSigErrorD.png");


    CreateHistogram(tree, "fGain", "(det==0)*(fSA_sig_err)", "hist", "Gain level vs sigma error in the Upstream sensor", "Noise level", "Sigma error", "fNoiseSigErrorU.png");
    CreateHistogram(tree, "fGain", "(det==1)*(fSA_sig_err)", "hist", "Gain level vs sigma error in the Downstream sensor", "Noise level", "Sigma error", "fNoiseSigErrorD.png");


    CreateHistogram(tree, "fOlScale", "(det==0)*(fSA_sig_err)", "hist", "OI Scale level vs sigma error in the Upstream sensor ", "Noise level", "Sigma error", "fNoiseSigErrorU.png");
    CreateHistogram(tree, "fOlScale", "(det==1)*(fSA_sig_err)", "hist", "OI Scale level vs sigma error in the Downstream sensor ", "Noise level", "Sigma error", "fNoiseSigErrorD.png");


    CreateHistogram(tree, "bunch", "(det==0)*fSA_mea", "hist", "Fit means per bunch in the Upstream sensor ", "Bunch ID", "Mean", "fNoiseSigErrorU.png");
    CreateHistogram(tree, "bunch", "(det==1)*fSA_mea", "hist", "Fit means per bunch in the Downstream sensor ", "Bunch ID", "Mean", "fNoiseSigErrorD.png");

    CreateHistogram(tree, "bunch", "(det==0)*fSA_sig", "hist", "Fit sigmas per bunch in the Upstream sensor ", "Bunch ID", "Sigma", "fNoiseSigErrorU.png");
    CreateHistogram(tree, "bunch", "(det==1)*fSA_sig", "hist", "Fit sigmas per bunch in the Downstream sensor ", "Bunch ID", "Sigma", "fNoiseSigErrorD.png");


    
    // Close the .root file
    file->Close();
    
    // Clean up
    delete file;
}

void OPTHistos() {
    CreateHistograms();
}