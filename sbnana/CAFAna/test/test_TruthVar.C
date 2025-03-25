#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

using namespace ana;

const TruthVar kTruthVar_NuE = SIMPLETRUTHVAR(E);

void test_TruthVar(){

  // Input
  SpectrumLoader loader("/pnfs/icarus/scratch/users/jskim/mc/NuMI_MC_Nu_Phase2/flatcaf/v09_72_00_03p01/230804_GENIESystSGV2_FixFSI/merged/flatcaf_0.root");

  // output; where plots are saved
  TString outputPath = "./test_TruthVar/";
  gSystem->mkdir(outputPath, kTRUE);

  // Binning
  const Binning binsEnergy = Binning::Simple(30, 0, 6);

  // All nu
  Spectrum *sTruthVar_NuE = new Spectrum("TruthVar_NuE", binsEnergy, loader, kTruthVar_NuE, kNoTruthCut, kNoSpillCut);

  // EnsembleSpectrum
  std::vector<TruthVar> vec_truthweights;
  for(int u=0; u<100; u++){
    std::string psetname = "GENIEReWeight_ICARUS_v2_multisim_ZExpAVariationResponse";
    vec_truthweights.push_back( GetTruthUniverseWeight(psetname, u) );
  }

  EnsembleSpectrum *es_TruthVar_Signal_NuE = new EnsembleSpectrum("ES_TruthVar_Signal_NuE", binsEnergy, loader, kTruthVar_NuE, kNoTruthCut, kNoSpillCut, vec_truthweights);

  loader.Go();

  TFile *f_out = new TFile(outputPath+"/output.root", "RECREATE");
  f_out->cd();

  TH1* h_TruthVar_NuE = sTruthVar_NuE->ToTH1(3e20);
  h_TruthVar_NuE->SetName("h_TruthVar_NuE");
  h_TruthVar_NuE->Write();

  TH1* h_TruthVar_NuE_FromES_Nominal = es_TruthVar_Signal_NuE->Nominal().ToTH1(3e20);
  h_TruthVar_NuE_FromES_Nominal->SetName("h_TruthVar_NuE_FromES_Nominal");
  h_TruthVar_NuE_FromES_Nominal->Write();
  TGraphAsymmErrors* gr_Error_FromES = es_TruthVar_Signal_NuE->ErrorBand(3e20);
  gr_Error_FromES->Write();
  f_out->Close();



}
