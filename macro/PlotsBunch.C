#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>
#include <TList.h>
#include <TObject.h>
#include <TStyle.h>
#include <iostream>

void PlotsBunch() {
       // Open the root file
    TFile *file = TFile::Open("BunchesOutput.root");
    if (!file) {
        std::cerr << "Failed to open the file" << std::endl;
        return;
    }

    // Array of canvas names
    const char* canvasNames[] = {
        "bunchIdVslgMeanUp", "bunchIdVshgMeanUp", "bunchIdVslgStdDevUp", "bunchIdVshgStdDevUp",
        "bunchIdVslgMeanDo", "bunchIdVshgMeanDo", "bunchIdVslgStdDevDo", "bunchIdVshgStdDevDo"
    };

    // Number of canvases
    const int numCanvases = sizeof(canvasNames) / sizeof(canvasNames[0]);

    for (int i = 0; i < numCanvases; i++) {
        // Retrieve the TCanvas
        TCanvas *c = (TCanvas*)file->Get(canvasNames[i]);
        if (!c) {
            std::cerr << "Failed to get the canvas: " << canvasNames[i] << std::endl;
            continue;
        }

        // Loop over the primitives in the TCanvas
        TIter next(c->GetListOfPrimitives());
        TObject *obj;
        while ((obj = next())) {
            if (obj->InheritsFrom(TGraph::Class())) {
                TGraph *g = (TGraph*)obj;
                g->SetMarkerStyle(24); // Use 24 for empty circles
                g->SetMarkerSize(1.0);
            }
        }

        // Update the TCanvas and save it as a .png file
        c->SetLeftMargin(0.15);
        c->Draw();
        c->SaveAs(("500kRunOutputs/" + std::string(canvasNames[i]) + ".png").c_str());    }

    file->Close();
}