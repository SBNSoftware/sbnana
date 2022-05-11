#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "cafanacore/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "TBox.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"

using namespace ana;

// Since I'm in the SBND group the SAM dataset works fine for me
const std::string sbnd_wildcard = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_flat_caf_sbnd";

void bar_sbnd()
{
  SpectrumLoader loader(sbnd_wildcard);

  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO

  const std::vector<std::pair<std::string, Cut>> cuts =
    {{"All Slices", kNoCut},
     {"NuScore", kSlcNuScoreCut},
     {"FV", kInFV},
     {"FlashMatchTime", kSlcFlashMatchTimeCut},
     {"FlashMatchScore", kSlcFlashMatchScoreCut},
     {"PrimaryMuonTrk", kHasPrimaryMuonTrk},
     {"CRTTrackAngle", kCRTTrackAngleCut},
     {"CRTHitDistance", kCRTHitDistanceCut}};

  const Cut kIsCosmic = SIMPLEVAR(truth.index) < 0;

  const HistAxis axCount("", Binning::Simple(1, 0, 1), Constant(.5));

  std::vector<Spectrum*> sNumuCC;
  std::vector<Spectrum*> sCosmic;
  std::vector<Spectrum*> sOther;

  Cut accumCut = kNoCut;
  for(auto& it: cuts){
    accumCut = accumCut && it.second;
    sNumuCC.push_back(new Spectrum(loader, axCount, kNumuSpillSel, accumCut && kIsNumuCC));
    sCosmic.push_back(new Spectrum(loader, axCount, kNumuSpillSel, accumCut && kIsCosmic));
    sOther.push_back(new Spectrum(loader, axCount, kNumuSpillSel, accumCut && !kIsNumuCC && !kIsCosmic));
  }

  loader.Go();

  const int Nbins = 4*cuts.size()+1;
  TH2* axes = new TH2F("", ";Neutrino candidates / 6.6#times10^{20} POT", 100, 3e4, 2e8, Nbins, 0, Nbins);
  axes->Draw();
  axes->GetXaxis()->CenterTitle();
  axes->GetYaxis()->SetLabelSize(0);
  gPad->SetLogx();
  gPad->SetLeftMargin(.3);

  TLegend* leg = new TLegend(.6, .15, .85, .35);
  leg->SetFillStyle(0);

  for(unsigned int i = 0; i < cuts.size(); ++i){
    const double nOther = sOther[i]->Integral(kPOTnominal);
    TBox* b = new TBox(0, Nbins-1-4*i, nOther, Nbins-2-4*i);
    b->SetFillColor(kGreen+2);
    b->Draw("lf same");
    if(i == 0) leg->AddEntry(b, "Other beam", "bf");
    TLatex* ltx = new TLatex(4e4, Nbins-1.5-4*i, TString::Format("%g", nOther));
    ltx->SetTextSize(.02);
    ltx->SetTextAlign(12);
    ltx->Draw();

    const double nCosmic = sCosmic[i]->Integral(kPOTnominal);
    b = new TBox(0, Nbins-2-4*i, nCosmic, Nbins-3-4*i);
    b->SetFillColor(kOrange+7);
    b->Draw("lf same");
    if(i == 0) leg->AddEntry(b, "Cosmic", "bf");
    ltx = new TLatex(4e4, Nbins-2.5-4*i, TString::Format("%g", nCosmic));
    ltx->SetTextSize(.02);
    ltx->SetTextAlign(12);
    ltx->Draw();

    const double nCC = sNumuCC[i]->Integral(kPOTnominal);
    b = new TBox(0, Nbins-3-4*i, nCC, Nbins-4-4*i);
    b->SetFillColor(kBlue);
    b->Draw("lf same");
    if(i == 0) leg->AddEntry(b, "#nu_{#mu} CC", "bf");
    ltx = new TLatex(4e4, Nbins-3.5-4*i, TString::Format("%g", nCC));
    ltx->SetTextSize(.02);
    ltx->SetTextAlign(12);
    ltx->SetTextColor(kWhite);
    ltx->Draw();

    TLatex* label = new TLatex(3e4, Nbins-2.5-4*i, (cuts[i].first+"  ").c_str());
    label->SetTextAlign(32);
    label->Draw();
  }

  leg->Draw();

  gPad->Print("bars_sbnd.png");
  gPad->Print("bars_sbnd.pdf");
}
