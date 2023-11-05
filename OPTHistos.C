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

// Function to create a graph from the TTree
void CreateGraph(TTree* tree, const char* var, const char* cut, 
                 const char* title, const char* xTitle, const char* yTitle, 
                 const char* fileName) {

    // Get the total number of entries in the tree
    Int_t nEntries = tree->Draw(var, cut, "goff");
    
    // Retrieve the arrays of x and y values from the tree
    Double_t *x = tree->GetV1();
    Double_t *y = tree->GetV2();
    
    // Create a new TGraph object using the arrays of x and y values
    TGraph *graph = new TGraph(nEntries, x, y);
    
    // Set the titles of the graph and its axes
    graph->SetTitle(title);
    graph->GetXaxis()->SetTitle(xTitle);
    graph->GetYaxis()->SetTitle(yTitle);
    
    // Create a new canvas to hold the graph
    TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 200, 10, 700, 500);
    
    // Draw the graph on the canvas
    graph->Draw("AP*");  // The option "AP" tells ROOT to draw axes and points
    
    // Save the canvas to a file
    canvas->SaveAs(fileName);
    
    // Clean up to prevent memory leaks
    delete graph;
    delete canvas;
}


void CreateHistograms() {
    // Open the .root file
    TFile* file = TFile::Open("dgtOut.root");
    
    // Get the TTree
    TTree* tree = (TTree*)file->Get("OPT");
    
    // Create histograms
    CreateHistogram(tree, "fNoise", "(det==0)", "hist", "Noise in the Upstream sensor Histogram", "Noise", "Counts", "fNoiseU.png");
    CreateHistogram(tree, "fNoise", "(det==1)", "hist", "Noise in the Downstream sensor Histogram", "Noise", "Counts", "fNoiseD.png");

    CreateGraph(tree, "fNoise:bunch", "(det==0)", "Noise per Bunch in the Upstream Histogram", "Bunch ID", "Noise", "fNoiseBunchU.png");
    CreateGraph(tree, "fNoise:bunch", "(det==1)", "Noise per Bunch in the Downstream Histogram", "Bunch ID", "Noise", "fNoiseBunchD.png");

    CreateHistogram(tree, "fOlScale", "(det==0)", "hist", "Fullscale Range in the Upstream Histogram", "Fullscale Range", "Counts", "fOlScaleU.png");
    CreateHistogram(tree, "fOlScale", "(det==1)", "hist", "Fullscale Range in the Downstream Histogram", "Fullscale Range", "Counts", "fOlScaleD.png");

    CreateGraph(tree, "fOlScale:bunch", "(det==0)", "Fullscale Range per Bunch in the Upstream Histogram", "Fullscale Range Bunch", "Counts", "fOlScaleBunchU.png");
    CreateGraph(tree, "fOlScale:bunch", "(det==1)", "Fullscale Range per Bunch in the Downstream Histogram", "Fullscale Range Bunch", "Counts", "fOlScaleBunchD.png");


    CreateHistogram(tree, "fGain", "(det==0)", "hist", "Gain in the Upstream Histogram", "Gain", "Counts", "fGainU.png");
    CreateHistogram(tree, "fGain", "(det==1)", "hist", "Gain in the Downstream Histogram", "Gain", "Counts", "fGainD.png");


    CreateGraph(tree, "fGain:bunch", "(det==0)", "Gain Bunch in the Upstream Histogram", "Gain Bunch", "Counts", "fGainBunchU.png");
    CreateGraph(tree, "fGain:bunch", "(det==1)", "Gain Bunch in the Downstream Histogram", "Gain Bunch", "Counts", "fGainBunchD.png");


    // Observables

    CreateGraph(tree, "(fSA_sig_err[0]):fNoise", "", "Noise level vs sigma error in the Upstream sensor ", "Noise level", "Sigma error", "fNoiseSigErrorU.png");
    CreateGraph(tree, "(fSA_sig_err[1]):fNoise", "", "Noise level vs sigma error in the Downstream sensor ", "Noise level", "Sigma error", "fNoiseSigErrorD.png");


    CreateGraph(tree, "(fSA_sig_err[0]):fGain", "", "Gain level vs sigma error in the Upstream sensor", "Noise level", "Sigma error", "fGainSigErrorU.png");
    CreateGraph(tree, "(fSA_sig_err[1]):fGain", "", "Gain level vs sigma error in the Downstream sensor", "Noise level", "Sigma error", "fGainSigErrorD.png");


    CreateGraph(tree, "(fSA_sig_err[0]):fOlScale", "", "OI Scale level vs sigma error in the Upstream sensor ", "Noise level", "Sigma error", "fOIScaleSigErrorU.png");
    CreateGraph(tree, "(fSA_sig_err[1]):fOlScale", "", "OI Scale level vs sigma error in the Downstream sensor ", "Noise level", "Sigma error", "fOIScaleSigErrorD.png");


    CreateHistogram(tree, "(fSA_mea[0]):bunch", "(det==0)", "hist", "Fit means per bunch in the Upstream sensor ", "Bunch ID", "Mean", "BunchMeanU.png");
    CreateHistogram(tree, "(fSA_mea[1]):bunch", "(det==1)", "hist", "Fit means per bunch in the Downstream sensor ", "Bunch ID", "Mean", "BunchMeanD.png");

    CreateHistogram(tree, "(fSA_sig[0]):bunch", "(det==0)", "hist", "Fit sigmas per bunch in the Upstream sensor ", "Bunch ID", "Sigma", "BunchSigU.png");
    CreateHistogram(tree, "(fSA_sig[1]):bunch", "(det==1)", "hist", "Fit sigmas per bunch in the Downstream sensor ", "Bunch ID", "Sigma", "BunchSigD.png");


    
    // Close the .root file
    file->Close();
    
    // Clean up
    delete file;
}

void OPTHistos() {
    CreateHistograms();
}