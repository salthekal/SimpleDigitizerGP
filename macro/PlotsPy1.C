#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TF1.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TF1.h>
#include <TLatex.h>

void ProcessHistogram(TFile *f, const char *hist_name, const char *fit_name, const char *output_file)
{
    // Check if the file is properly opened
    if (!f || f->IsZombie()) {
        std::cout << "Error: File not opened successfully" << std::endl;
        return;
    }
    
    // Get histogram
    TH1D *h = (TH1D *)f->Get(hist_name);
    if (!h) {
        std::cout << "Error: Histogram " << hist_name << " not found in file" << std::endl;
        return;
    }
    
    // Check if the histogram has any entries
    if (h->GetEntries() == 0) {
        std::cout << "Warning: Histogram " << hist_name << " is empty" << std::endl;
        // Optionally, you can return here if an empty histogram should not be processed further
        // return;
    }
    
    // Create a canvas
    TCanvas *c = new TCanvas();
    gStyle->SetOptStat(1001111);

    // Draw histogram with errors
    h->Draw("E0");

    // Check histogram range
    double histMin = h->GetXaxis()->GetXmin();
    double histMax = h->GetXaxis()->GetXmax();
    if (histMin > 60 || histMax < 160)
    {
        std::cout << "Error: Histogram range does not include fit range" << std::endl;
        return;
    }

    // Check data in fit range
    double integral = h->Integral(h->FindBin(60), h->FindBin(160));
    if (integral == 0)
    {
        std::cout << "Warning: No data in fit range" << std::endl;
        // Optionally, return here if no data in fit range should not be processed further
        // return;
    }
    // Fit the graph
    //TF1 *fit = new TF1("fit", fit_name, histMin, histMax);
   // Fit the histogram with a predefined exponential function
    h->Fit("gaus");

    // Check the fit status
    Int_t fitResult = h->Fit("gaus", "S");
    if (fitResult != 0)
    {
        std::cout << "Error: Fit did not converge" << std::endl;
        return;
    }

    // Draw the histogram and the fit
    h->Draw("SAME");
    h->GetFunction("gaus")->SetLineColor(kRed); // Optional: change the color of the fit line
    h->GetFunction("gaus")->Draw("same");       // Enable the display of fit parameters in the stats box
    gStyle->SetOptFit(1111);
    c->Update();

    // Adjust stats box position
    TPaveStats *st = (TPaveStats *)c->GetPrimitive("stats");
    if (st != nullptr)
    {
        st->SetX1NDC(0.71);
        st->SetX2NDC(0.9);
        st->SetY1NDC(0.71);
        st->SetY2NDC(0.9);
        st->Draw();
    }
    else
    {
        std::cout << "Error: Stats object not found" << std::endl;
    }

    c->Draw();
    c->SaveAs(output_file);
}

void ProcessHistogramWithoutFit(TFile *f, const char *hist_name, const char *output_file, bool special_margin = false)
{
    TCanvas *c = new TCanvas();
    if (special_margin) {
        c->SetLeftMargin(0.15);
    }

    TH1D *h = (TH1D *)f->Get(hist_name);
    h->Draw("hist");
    gStyle->SetOptStat(1001111);

    TPaveStats *st = (TPaveStats *)h->FindObject("stats");
    st->SetX1NDC(0.71);
    st->SetX2NDC(0.9);
    st->SetY1NDC(0.71);
    st->SetY2NDC(0.9);
    st->Draw();

    c->SaveAs(output_file);
}

void PlotsPy1()
{

    TFile *f = TFile::Open("digitizationGamma1M.root"); // replace with your file

    ProcessHistogramWithoutFit(f, "hist_chgDepProfileY;1", "10MGammaRun/chgDepUp.png");
    ProcessHistogramWithoutFit(f, "hist_chgDepProfileX;2", "10MGammaRun/chgDepDo.png");
    ProcessHistogramWithoutFit(f, "hist_chgProjProfileY;1", "10MGammaRun/chgProjUp.png");
    ProcessHistogramWithoutFit(f, "hist_chgProjProfileX;2", "10MGammaRun/chgProjDo.png");
    ProcessHistogramWithoutFit(f, "hist_chgProjFEProfileX;2", "10MGammaRun/chgProjFEDo.png", true);
    ProcessHistogramWithoutFit(f, "hist_chgProjFEProfileY;1", "10MGammaRun/chgProjFEUp.png", true);
    ProcessHistogramWithoutFit(f, "hist_chgDepProfileY_strip;1", "10MGammaRun/chgDepStripUp.png");
    ProcessHistogramWithoutFit(f, "hist_chgDepProfileX_strip;2", "10MGammaRun/chgDepStripDo.png");
    ProcessHistogramWithoutFit(f, "hist_chgProjProfileY_strip;1", "10MGammaRun/chgProjStripUp.png");
    ProcessHistogramWithoutFit(f, "hist_chgProjProfileX_strip;2", "10MGammaRun/chgProjStripDo.png");
    ProcessHistogramWithoutFit(f, "hist_lgADCinVoltProfX_strip;2", "10MGammaRun/lgADCinVoltDo.png", true);
    ProcessHistogramWithoutFit(f, "hist_lgADCinVoltProfY_strip;1", "10MGammaRun/lgADCinVoltUp.png", true);
    ProcessHistogramWithoutFit(f, "hist_hgADCinVoltProfX_strip;2", "10MGammaRun/hgADCinVoltDo.png", true);
    ProcessHistogramWithoutFit(f, "hist_hgADCinVoltProfY_strip;1", "10MGammaRun/hgADCinVoltUp.png", true);
    ProcessHistogramWithoutFit(f, "hist_chgProjFEProfileX_strip;2", "10MGammaRun/chgProjFEStripDo.png", true);
    ProcessHistogramWithoutFit(f, "hist_chgProjFEProfileY_strip;1", "10MGammaRun/chgProjFEStripUp.png", true);
    ProcessHistogramWithoutFit(f, "hist_lgADCctsProfX_strip;2", "10MGammaRun/lgADCctsProfDo.png", true);
    ProcessHistogramWithoutFit(f, "hist_lgADCctsProfY_strip;1", "10MGammaRun/lgADCctsProfUp.png", true);
    ProcessHistogramWithoutFit(f, "hist_hgADCctsProfX_strip;2", "10MGammaRun/hgADCctsProfDo.png", true);
    ProcessHistogramWithoutFit(f, "hist_hgADCctsProfY_strip;1", "10MGammaRun/hgADCctsProfUp.png", true);
    ProcessHistogramWithoutFit(f, "hist_chgProjProfileY_strip;3", "10MGammaRun/chgProjCrossProfUp.png");
    ProcessHistogramWithoutFit(f, "hist_chgProjProfileX_strip;3", "10MGammaRun/chgProjCrossProfDo.png");
    ProcessHistogram(f, "hist_lgADCctsProfX_strip;2", "pow", "10MGammaRun/lgADCctsProfDoFit.png");
    ProcessHistogram(f, "hist_lgADCctsProfY_strip;1", "pow", "10MGammaRun/lgADCctsProfUpFit.png");
    ProcessHistogram(f, "hist_hgADCctsProfY_strip;1", "pow", "10MGammaRun/hgADCctsProfUpFit.png");
    ProcessHistogram(f, "hist_hgADCctsProfX_strip;2", "pow", "10MGammaRun/hgADCctsProfDoFit.png");
    
    //////////////////////////////////////////////////////////////////////////////////////////////

    // Create a table for the resutls of the fits of the ADC profiles for each bunch
    // Create a canvas for the table
    TCanvas *canvas = new TCanvas("canvas", "Fit Parameters", 800, 600);

    // Clear the canvas
    canvas->Clear();

    // Create a TPaveText object to display the table
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    pt->SetTextAlign(12);     // Align the text to the left
    pt->SetFillColor(kWhite); // Set the background color to white

    // @TODO create histograms by storing the error of gaus fit of ADC profiles and plotting it as a function of FE noise sigma (x-axis)

    
    // Create an ofstream object for writing to a file
    std::ofstream outFile("FitResults.txt");

    // Generate histogram names and fit them
    for (int i = 1; i <= 10; ++i) { // assuming 10 histograms
        for (int j = 1; j <= 4; ++j) { // assuming 4 versions per histogram
            std::string histName = "hist_" + std::to_string(i) + ";" + std::to_string(j);
        // Get the histogram
        TH1F *hist = (TH1F *)f->Get(histName.c_str());

        if (hist)
        {
         // Fit the histogram with a Gaussian
        hist->Fit("gaus", "SAME"); // "Q" for quiet mode

        // Get the fit parameters
        double constant = fit->GetParameter(0);
        double mean = fit->GetParameter(1);
        double stdDev = fit->GetParameter(2);

        // Get the errors of the fit parameters
        double constantError = fit->GetParError(0);
        double meanError = fit->GetParError(1);
        double stdDevError = fit->GetParError(2);

        // Get the fit quality (chi-square / number of degrees of freedom)
        double chi2NDF = fit->GetChisquare() / fit->GetNDF();

        // Add the parameters to the table
        TString line = TString::Format("Histogram %-10s : Constant = %.2f ± %.2f, Mean = %.2f ± %.2f, Std Dev = %.2f ± %.2f, χ²/NDF = %.2f",
                                       histName.c_str(), constant, constantError, mean, meanError, stdDev, stdDevError, chi2NDF);
        pt->AddText(line);

         // Write stdDev and stdDevError to the file
        outFile << "Histogram " << histName << ": Std Dev = " << stdDev
             << ", Std Dev Error = " << stdDevError << std::endl;
        }
        else
        {
            std::cerr << "Failed to retrieve histogram: " << histName << std::endl;
        }
    } 
    // Close the output file
    outFile.close();
    
    // Draw the table
    pt->Draw();

    // Update the canvas
    canvas->Update();

    // Save the canvas as a .png file
    canvas->SaveAs("10MGammaRun/FitParameters.png");
    f->Close();
////////////////////////////////////////////////////////////////////////////

}