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

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"

using namespace ana;

// Since I'm in the SBND group the SAM dataset works fine for me
const std::string sbnd_wildcard = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_flat_caf_sbnd";

const std::string state_fname = "test_ensemble_state.root";

void test_ensemble(bool reload = false)
{
  if(reload || TFile(state_fname.c_str()).IsZombie()){
    SpectrumLoader loader(sbnd_wildcard);

    const Var kTrueE = SIMPLEVAR(truth.E);

    const Cut kNumuSel = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut;

    const Binning binsEnergy = Binning::Simple(30, 0, 3);

    const HistAxis axEnergy("True energy (GeV)", binsEnergy, kTrueE);

    std::vector<Var> weis;
    weis.reserve(100);

    // We need these essentially as the list of knob names
    const std::vector<const ISyst*>& systs = GetSBNGenieWeightSysts();
    for(int i = 0; i < 100; ++i) weis.push_back(GetUniverseWeight(systs, i));

    EnsembleSpectrum sCC(loader, axEnergy, kNoSpillCut, kNumuSel && kIsNumuCC, weis);
    EnsembleSpectrum sNC(loader, axEnergy, kNoSpillCut, kNumuSel && kIsNC, weis);

    loader.Go();

    TFile fout(state_fname.c_str(), "RECREATE");
    sCC.SaveTo(fout.mkdir("cc"));
    sNC.SaveTo(fout.mkdir("nc"));
  }

  TFile fin(state_fname.c_str());
  EnsembleSpectrum* sCC = LoadFrom<EnsembleSpectrum>(fin.GetDirectory("cc")).release();
  EnsembleSpectrum* sNC = LoadFrom<EnsembleSpectrum>(fin.GetDirectory("nc")).release();

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
