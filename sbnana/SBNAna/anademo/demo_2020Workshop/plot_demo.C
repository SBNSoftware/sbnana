#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"

#include "demo.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"


using namespace ana;

void plot_demo(std::string inFile)
{
  // Arbitrary POT to scale all plots to
  // Spill tree currently not being filled, so there will be a warning
  double POT = 1E20;

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();

  // I want to make a plot for each var
  for(unsigned int iVar = 0; iVar < kNVar; ++iVar){
    
    TCanvas *c = new TCanvas(plots[iVar].suffix.c_str(),plots[iVar].suffix.c_str());

    std::vector<TH1*> hists;
    TLegend *l = new TLegend(0.15, 0.6, 0.4, 0.8);
    l->SetFillStyle(0);

    for(unsigned int jSel = 0; jSel < kNSel; ++jSel){
      std::string mysuffix = sels[jSel].suffix + "_" + plots[iVar].suffix;

      Spectrum *spec = LoadFromFile<Spectrum>(inFile, mysuffix).release();
      TH1* h = spec->ToTH1(POT);
      h->SetLineColor(sels[jSel].color);

      hists.push_back(h);
      l->AddEntry(h, sels[jSel].label.c_str(), "l");
    } // iSel

    hists[0]->GetXaxis()->CenterTitle();
    hists[0]->GetYaxis()->CenterTitle();

    hists[0]->Draw("hist");
    for(unsigned int jSel = 0; jSel < kNSel; ++jSel) 
      hists[jSel]->Draw("hist same");
    hists[0]->Draw("hist same");
    l->Draw();

    c->Print(("plots/"+plots[iVar].suffix+".pdf").c_str());

  } // iVar
}
