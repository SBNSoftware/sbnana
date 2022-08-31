// Recommend running with a large --stride argument to just use the first file
// in the dataset

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"

using namespace ana;

const std::string wildcard = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_37_01_03p02/MultiSigmaAdded_RenamedPSet/*.root";

const std::string state_fname = "test_ensemble_state.root";

void test_ensemble(bool reload = false)
{
  if(reload || TFile(state_fname.c_str()).IsZombie()){
    SpectrumLoader loader(wildcard);

    const Var kTrueE = SIMPLEVAR(truth.E);

    const Cut kNumuSel = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut;

    const Binning binsEnergy = Binning::Simple(30, 0, 3);

    const HistAxis axEnergy("True energy (GeV)", binsEnergy, kTrueE);

    std::vector<Weight> weis;
    weis.reserve(101);

    weis.push_back(kUnweighted); // nominal
    for(int i = 0; i < 99; ++i) weis.push_back(GetUniverseWeight("multisim_Genie", i));

    EnsembleSpectrum sCC(loader.Slices().Ensemble(weis)[kNumuSel][kIsNumuCC], axEnergy);
    EnsembleSpectrum sNC(loader.Slices().Ensemble(weis)[kNumuSel][kIsNC],     axEnergy);

    loader.Go();

    TFile fout(state_fname.c_str(), "RECREATE");
    sCC.SaveTo(&fout, "cc");
    sNC.SaveTo(&fout, "nc");
  }

  TFile fin(state_fname.c_str());
  EnsembleSpectrum* sCC = LoadFrom<EnsembleSpectrum>(&fin, "cc").release();
  EnsembleSpectrum* sNC = LoadFrom<EnsembleSpectrum>(&fin, "nc").release();

  sCC->Nominal().ToTH1(kPOTnominal, kRed)->Draw("hist");
  sNC->Nominal().ToTH1(kPOTnominal, kBlue)->Draw("hist same");

  for(unsigned int i = 0; i < sCC->NUniverses(); ++i){
    sCC->Universe(i).ToTH1(kPOTnominal, kRed-10)->Draw("hist same");
    sNC->Universe(i).ToTH1(kPOTnominal, kBlue-10)->Draw("hist same");
  }

  // Redraw the nominals over the top
  sCC->Nominal().ToTH1(kPOTnominal, kRed)->Draw("hist same");
  sNC->Nominal().ToTH1(kPOTnominal, kBlue)->Draw("hist same");

  gPad->Print("test_ensemble.pdf");

  new TCanvas;

  TGraphAsymmErrors* bandCC = sCC->ErrorBand(kPOTnominal);
  DrawErrorBand(sCC->Nominal().ToTH1(kPOTnominal, kRed), bandCC);
  TGraphAsymmErrors* bandNC = sNC->ErrorBand(kPOTnominal);
  DrawErrorBand(sNC->Nominal().ToTH1(kPOTnominal, kBlue), bandNC, kBlue, .5);

  gPad->Print("test_ensemble_band.pdf");

  new TCanvas;

  EnsembleRatio ratio = *sCC / *sNC;
  TH1* hratio = ratio.Nominal().ToTH1();
  hratio->Draw("hist ][");
  hratio->GetYaxis()->SetTitle("CC / NC ratio");
  hratio->GetYaxis()->SetRangeUser(0, 10);
  for(unsigned int i = 0; i < ratio.NUniverses(); ++i){
    ratio.Universe(i).ToTH1(kGray)->Draw("hist ][ same");
  }
  hratio->Draw("hist ][ same");

  gPad->Print("test_ensemble_ratio.pdf");

  new TCanvas;
  TGraphAsymmErrors* bandRatio = ratio.ErrorBand();
  DrawErrorBand(hratio, bandRatio, kGray);

  gPad->Print("test_ensemble_ratio_band.pdf");
}
