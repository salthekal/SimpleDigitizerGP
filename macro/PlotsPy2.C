#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TObject.h"


// Helper function to draw and save histograms
void draw_and_save(TH1* hist, TObject* stats, const char* output_name) {
    hist->Draw("hist");
    gStyle->SetOptStat(1001111);
    hist->SetLeftMargin(0.15);
    TPaveStats *st = (TPaveStats*)stats->FindObject("stats");
    st->SetX1NDC(0.71);
    st->SetX2NDC(0.9);
    st->SetY1NDC(0.71);
    st->SetY2NDC(0.9);
    st->Draw();
    hist->SaveAs(output_name);
}

void PlotsPy2() {
    TFile* f = TFile::Open("BunchesOutput.root");

    // Draw and save the histograms
    draw_and_save(lMeanHistUp, lgadc_meansUp, "1MGammaRun/lMeanHistUp.png");
    draw_and_save(hMeanHistUp, hgadc_meansUp, "1MGammaRun/hMeanHistUp.png");
    draw_and_save(lStdDevsHistUp, lgadc_StdDevsUp, "1MGammaRun/lStdDevsHistUp.png");
    draw_and_save(hStdDevsHistUp, hgadc_StdDevsUp, "1MGammaRun/hStdDevsHistUp.png");
    draw_and_save(lMeanHistDo, lgadc_meansDo, "1MGammaRun/lMeanHistDo.png");
    draw_and_save(hMeanHistDo, hgadc_meansDo, "1MGammaRun/hMeanHistDo.png");
    draw_and_save(lStdDevsHistDo, lgadc_StdDevsDo, "1MGammaRun/lStdDevsHistDo.png");
    draw_and_save(hStdDevsHistDo, hgadc_StdDevsDo, "1MGammaRun/hStdDevsHistDo.png");

    // For these, we don't have a stats object associated, so we simply draw and save.
    lgADCVoltVsADCUp->Draw("hist");
    lgADCVoltVsADCUp->SaveAs("1MGammaRun/lgADCVoltVsADCUp.png");
    hgADCVoltVsADCUp->Draw("hist");
    hgADCVoltVsADCUp->SaveAs("1MGammaRun/hgADCVoltVsADCUp.png");
    lgADCVoltVsADCDo->Draw("hist");
    lgADCVoltVsADCDo->SaveAs("1MGammaRun/lgADCVoltVsADCDo.png");
    hgADCVoltVsADCDo->Draw("hist");
    hgADCVoltVsADCDo->SaveAs("1MGammaRun/hgADCVoltVsADCDo.png");

}