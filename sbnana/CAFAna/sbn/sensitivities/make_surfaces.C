#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Vars/FitVarsSterileApprox.h"
#include "sbnana/CAFAna/Prediction/PredictionInterp.h"
#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"
#include "sbnana/CAFAna/Experiment/MultiExperimentSBN.h"
#include "sbnana/CAFAna/Experiment/CountingExperiment.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/Analysis/Surface.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/BoosterFluxSysts.h"
#include "sbnana/CAFAna/Systs/EnergySysts.h"
#include "sbnana/CAFAna/Systs/Systs.h"
//#include "sbnana/CAFAna/Systs/Systs.h"
#include "sbnana/CAFAna/Systs/UniverseOracle.h"
using namespace ana;

#include "OscLib/IOscCalc.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"

#include <vector>
#include <fstream>

// This file includes the detector systematics in my local directory.
// When I commit those to the repository this line will need to change.
// #include "mySysts.h"

//const double sbndPOT = kPOTnominal;
//const double icarusPOT = kPOTnominal;
const double sbndPOT = 10e20;
const double icarusPOT = 2.5e20;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void make_surfaces(const std::string anatype = numuStr)
{
  const auto xsec_systs = GetSBNGenieWeightSysts();
  const auto beam_systs = GetSBNBoosterFluxWeightSysts();
  const auto flux_had_systs = GetBoosterFluxHadronSysts(30);

  std::vector<const ISyst*> flux_systs = beam_systs;
  flux_systs.insert(std::end(flux_systs), std::begin(flux_had_systs), std::end(flux_had_systs));
  
  std::vector<const ISyst*> xsec_flux = xsec_systs;
  xsec_flux.insert(std::end(xsec_flux), std::begin(flux_systs), std::end(flux_systs));

  const char* name_in;
  const char* name_out;
  if (anatype == numuStr) {
    name_in = "cafe_state_syst_numu.root";
    name_out = "surfaces_numu.root";
  }
  else if (anatype == nueStr) {
    name_in = "cafe_state_syst_nue.root";
    name_out = "surfaces_nue.root";
  }
  else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }

  TFile fin(name_in);
  TFile fout(name_out,"RECREATE");

  //PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();

  //Cosmics need to be added in separately 
  //Spectrum* cosmics_np = LoadFrom<Spectrum>(fin.GetDirectory("cosmics_np_spec")).release();
  //Spectrum* intime_np = LoadFrom<Spectrum>(fin.GetDirectory("intime_np_spec")).release();

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();
  OscCalcSterileApproxAdjustable* seed = DefaultSterileApproxCalc();

  if (anatype == nueStr) {
    seed->calc.SetSinSq2ThetaMuE(1e-3);
    seed->calc.SetDmsq(1.32);
    //p_nd->SetOscSeed(seed);
    p_fd->SetOscSeed(seed);
  }

  const IFitVar* SterileTh;
  double thlim = 1e-3;
  if (anatype == numuStr) {
    SterileTh = &kFitSinSq2ThetaMuMu;
  }
  else {
    SterileTh = &kFitSinSq2ThetaMuE;
    thlim = 5e-5;
  }
  const FitAxis kAxForTh(SterileTh, 40, thlim, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 40, 1e-2, 1e2, true);

  osc::NoOscillations noosc;
  //Spectrum data_nd = p_nd->Predict(&noosc).FakeData(sbndPOT);
  Spectrum data_fd = p_fd->Predict(&noosc).FakeData(2.5e20);
  // With oot cosmics
  //Spectrum data_fd_np = (p_fd_np->Predict(&noosc)+*cosmics_np).FakeData(2.5e20);
  
  double targetLivetimeFD = 2.5e20 / (5e12);

  // Add in in-time cosmics. This OverideLivetime call is needed to force
  // the POT/Livetime scaling to do the right thing when we add in in-time cosmics
  //data_fd_np.OverrideLivetime(targetLivetimeFD - data_fd_np.Livetime());
  //data_fd_np += *intime_np;
  //data_fd_np.OverrideLivetime(targetLivetimeFD);

  //SingleSampleExperiment expt_nd(p_nd, data_nd);
  SingleSampleExperiment expt_fd(p_fd, data_fd);
  // Pass cosmics spectra to SingleSampleExperiment
  //SingleSampleExperiment expt_fd_np(p_fd_np, data_fd_np, *intime_np, *cosmics_np);

  //MultiExperimentSBN fd_nd({&expt_nd, &expt_fd}, {kSBND, kICARUS});

  //p_nd->MinimizeMemory();
  p_fd->MinimizeMemory();
  //p_fd_np->MinimizeMemory();
  
  //Surface surf_nd_fd(&fd_nd, calc, kAxForTh, kAxDmSq);
  //Surface surf_nom_nd(&expt_nd, calc, kAxForTh, kAxDmSq);
  Surface surf_nom_fd(&expt_fd, calc, kAxForTh, kAxDmSq);
  //Surface surf_nom_fd_np(&expt_fd_np, calc, kAxForTh, kAxDmSq);

  fout.mkdir("exclusion");
  fout.cd("exclusion");    

  //surf_nom_nd.SaveTo(gDirectory->mkdir("nom_nd"));
  surf_nom_fd.SaveTo(gDirectory->mkdir("nom_fd"));
  //surf_nom_fd_np.SaveTo(gDirectory->mkdir("nom_fd_np"));
  //surf_nd_fd.SaveTo(gDirectory->mkdir("nom_nd_fd"));

  std::map<std::string, std::vector<const ISyst*>> slists;
  slists["xsec"] = xsec_systs;
  slists["flux"] = flux_systs;
  slists["xsec_flux"] = xsec_flux;

  for(const auto& [suffix, slist]: slists){
    //Surface surf_syst_nd_fd(&fd_nd, calc, kAxForTh, kAxDmSq, {}, slist);
    //Surface surf_syst_nd(&expt_nd, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_fd(&expt_fd, calc, kAxForTh, kAxDmSq, {}, slist);
    
    //surf_syst_nd_fd.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));
    //surf_syst_nd.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
    surf_syst_fd.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
  }

  // Allowed Region

//  OscCalcSterileApproxAdjustable* calc2 = DefaultSterileApproxCalc();
//  if (anatype == numuStr) {
//    calc2->calc.SetSinSq2ThetaMuMu(4*0.135*0.135*(1-0.135*0.135));
//    calc2->calc.SetDmsq(1.32);
//  }
//  else {
//    calc2->calc.SetSinSq2ThetaMuMu(4*0.135*0.135*(1-0.135*0.135));
//    calc2->calc.SetDmsq(1.32);
//    std::cout << "WARNING WRONG INJECTED VALUES FIX ME!!!!" << std::endl;
//  }
//
//  const Spectrum data_nd2 = p_nd->Predict(calc2).FakeData(sbndPOT);
//  const Spectrum data_fd2 = p_fd->Predict(calc2).FakeData(icarusPOT);
//
//  SingleSampleExperiment expt_nd2(p_nd, data_nd2);
//  SingleSampleExperiment expt_fd2(p_fd, data_fd2);
//
//  MultiExperimentSBN fd_nd2({&expt_nd2, &expt_fd2}, {kSBND, kICARUS});
//
//  //Surface surf_nd_fd2(&fd_nd2, calc2, kAxForTh, kAxDmSq);
//  //Surface surf_nom_nd2(&expt_nd2, calc2, kAxForTh, kAxDmSq);
//  //Surface surf_nom_fd2(&expt_fd2, calc2, kAxForTh, kAxDmSq);
//    
//  //fout.cd("..");
//  //fout.mkdir("allowed");
//  //fout.cd("allowed");
//
//  //surf_nom_nd2.SaveTo(gDirectory->mkdir("nom_nd"));
//  //surf_nom_fd2.SaveTo(gDirectory->mkdir("nom_fd"));
//  //surf_nd_fd2.SaveTo(gDirectory->mkdir("nom_nd_fd"));
//
//  for(const std::vector<const ISyst*> slist: slists){
//    //Surface surf_syst_nd_fd2(&fd_nd2, calc2, kAxForTh, kAxDmSq, {}, slist);
//    //Surface surf_syst_nd2(&expt_nd2, calc2, kAxForTh, kAxDmSq, {}, slist);
//    //Surface surf_syst_fd2(&expt_fd2, calc2, kAxForTh, kAxDmSq, {}, slist);
//
//    std::string suffix = "prop_systs";
//
//    //surf_syst_nd2.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
//    //surf_syst_fd2.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
//    //surf_syst2.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
//    //surf_syst_nd_fd2.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));
//
//  } // end for s
}
