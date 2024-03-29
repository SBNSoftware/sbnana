// Get Numbers for various interaction types and flavours
// cafe event_numbers.C

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

// #include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/Analysis/Calcs.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "OscLib/OscCalcSterile.h"
#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"

using namespace ana;

// Put an "SBND simulation" tag in the corner                                                       
void Experiment(std::string expt)
{
  TLatex* prelim = new TLatex(.9, .95, (expt+" Simulation").c_str());
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

void event_numbers_nue(const std::string expt = "SBND")
{
  const std::string dir = "/sbnd/data/users/bzamoran/Modern_6thNov2019/";
  const std::string fnameBNB  = dir + "output_SBNOsc_NueSelection_Modern_" + expt + "_Numu.flat.root";
  const std::string fnameInt  = dir + "output_SBNOsc_NueSelection_Modern_" + expt + "_Int.flat.root";
  const std::string fnameSwap = dir + "output_SBNOsc_NueSelection_Modern_" + expt + "_Osc.flat.root";

  // Source of events
  Loaders loaders;

  loaders.SetLoaderPath( fnameBNB,   Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  loaders.SetLoaderPath( fnameInt,   Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
  loaders.SetLoaderPath( fnameSwap,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

  // Vars and cuts

  // a struct to hold the different flavours and signs
  const unsigned int kNumFlavours = 3;
  const unsigned int kNumSigns    = 2;
  const unsigned int kNumCurrents = 2;

  struct Flavour
  {
    ana::Flavors::Flavors_t flav;
    string label;
  };
  Flavour flav[kNumFlavours] = {
    {ana::Flavors::kAllNuMu, "numu"},
    {ana::Flavors::kAllNuE,  "nue"},
    {ana::Flavors::kAll,  "all"}
  };
  struct Signtype
  {
    ana::Sign::Sign_t sign;
    string label;
  };
    Signtype sign[kNumSigns] = {
    {ana::Sign::kNu,     ""},
    {ana::Sign::kAntiNu, "anti"}
  };
  struct Currtype
  {
    ana::Current::Current_t curr;
    string label;
  };
    Currtype curr[kNumCurrents] = {
    {ana::Current::kCC, "CC"},
    {ana::Current::kNC, "NC"}
  };

  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  const Var kWeight = SIMPLEVAR(reco.weight);

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  const Var kCC = SIMPLEVAR(truth[0].neutrino.iscc);
  const Cut kIsCC = kCC > 0;

  // LSND best fit
  OscCalcSterileApproxAdjustable* osc_lsnd = DefaultSterileApproxCalc();
  osc_lsnd->SetL(expt == "SBND" ? kBaselineSBND : kBaselineIcarus);
  osc_lsnd->calc.SetSinSq2ThetaMuE(0.003);
  osc_lsnd->calc.SetDmsq(1.2);

  // The main predictions
  PredictionNoExtrap* pred = new PredictionNoExtrap(loaders, axEnergy, kNoCut, kNoShift, kWeight);

  // GO!
  loaders.Go();

  const double pot = kPOTnominal;

  TH1* h_numu = pred->PredictComponent(osc_lsnd, Flavors::kAllNuMu,   Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* h_NC   = pred->PredictComponent(osc_lsnd, Flavors::kAll,       Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* h_int  = pred->PredictComponent(osc_lsnd, Flavors::kNuEToNuE,  Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* h_osc  = pred->PredictComponent(osc_lsnd, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);

  std::cout << "Numu: " << h_numu->Integral() << std::endl;
  std::cout << "NC:   " << h_NC->Integral() << std::endl;
  std::cout << "Int:  " << h_int->Integral() << std::endl;
  std::cout << "Osc:  " << h_osc->Integral() << std::endl;

  TFile* fOutput = new TFile(("output/output_nue_"+expt+".root").c_str(),"RECREATE");

  TCanvas* c1 = new TCanvas("c1","c1");

  h_numu->SetLineColor(kGreen+2);
  h_numu->SetMarkerColor(kGreen+2);
  h_numu->Write((expt+"_numu").c_str());   

  h_NC->SetLineColor(kBlack);
  h_NC->SetMarkerColor(kBlack);
  h_NC->Write((expt+"_NC").c_str()); 

  h_int->SetLineColor(kRed);
  h_int->SetMarkerColor(kRed);
  h_int->Write((expt+"_int").c_str()); 

  h_osc->SetLineColor(kBlue);
  h_osc->SetMarkerColor(kBlue);
  h_osc->Write((expt+"_osc").c_str()); 

  fOutput->Close();

}
