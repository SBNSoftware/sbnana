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

const SpillMultiVar kNuEs([](const caf::SRSpillProxy *sr) -> vector<double> {
  std::vector<double> rets;
  for(const auto& nu: sr->mc.nu){
    rets.push_back(nu.E);
  }
  return rets;
});

const TruthVar kTruthVar_NuE = SIMPLETRUTHVAR(E);

void test_TruthVar(){

  // Input
  SpectrumLoader loader("/pnfs/icarus/scratch/users/jskim/mc/NuMI_MC_Nu_Phase2/flatcaf/v09_72_00_03p01/230919_Reproc_SkipCorrupt/flatcaf_0.root");

  // output; where plots are saved
  TString outputPath = "./test_TruthVar/";
  gSystem->mkdir(outputPath, kTRUE);

  // Binning
  const Binning binsEnergy = Binning::Simple(30, 0, 6);

  // Using SpillMultiVar
  Spectrum *sNuE = new Spectrum("NuE", binsEnergy, loader, kNuEs, kNoSpillCut);

  // Using TruthVar
  Spectrum *sTruthVar_NuE = new Spectrum("TruthVar_NuE", binsEnergy, loader, kTruthVar_NuE, kNoTruthCut, kNoSpillCut);

  // +1 sigma RPA_CCQE dial
  ISyst* iS = new SBNWeightSyst("GENIEReWeight_ICARUS_v2_multisigma_RPA_CCQE");
  Spectrum *sTruthVar_NuE_1up = new Spectrum("TruthVar_NuE_1up", binsEnergy, loader, kTruthVar_NuE, kNoTruthCut, kNoSpillCut, SystShifts(iS,+1));

  loader.Go();

  TFile *f_out = new TFile(outputPath+"/output.root", "RECREATE");
  f_out->cd();

  TH1* h_NuE = sNuE->ToTH1(3e20);
  h_NuE->SetName("h_NuE");
  h_NuE->Write();

  TH1* h_TruthVar_NuE = sTruthVar_NuE->ToTH1(3e20);
  h_TruthVar_NuE->SetName("h_TruthVar_NuE");
  h_TruthVar_NuE->Write();

  TH1* h_TruthVar_NuE_1up = sTruthVar_NuE_1up->ToTH1(3e20);
  h_TruthVar_NuE_1up->SetName("h_TruthVar_NuE_1up");
  h_TruthVar_NuE_1up->Write();




}
