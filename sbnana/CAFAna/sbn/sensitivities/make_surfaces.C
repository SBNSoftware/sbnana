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

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void make_surfaces(const std::string anatype = numuStr)
{
  const auto genie_systs = GetSBNGenieWeightSysts();
  const auto flux_systs = GetSBNFluxHadronSysts(30);
  const auto energy_systs = GetEnergySysts();
  
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

  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
  PredictionInterp* p_ub = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_ub")).release();

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();
  OscCalcSterileApproxAdjustable* seed = DefaultSterileApproxCalc();

  //JL - try different values here, how much does it matter what values we choose, does this
  //have any impact on disappearance? 
  if (anatype == nueStr) {
    seed->calc.SetSinSq2ThetaMuE(1e-3);
    seed->calc.SetDmsq(1.32);
    p_nd->SetOscSeed(seed);
    p_fd->SetOscSeed(seed);
    p_ub->SetOscSeed(seed);
  }
    
  //Define fit axes
  //smaller fit axes for no
  const IFitVar* SterileTh;
  double thlim = 1e-3;
  if (anatype == numuStr) {
    SterileTh = &kFitSinSq2ThetaMuMu;
  }
  else {
    SterileTh = &kFitSinSq2ThetaMuE;
    thlim = 5e-5;
  }
  const FitAxis kAxForTh(SterileTh, 10, thlim, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 10, 1e-2, 1e2, true);

  // We'll call zero nominal
  const Spectrum data_nd = p_nd->Predict(calc).FakeData(sbndPOT);
  const Spectrum data_fd = p_fd->Predict(calc).FakeData(icarusPOT);
  const Spectrum data_ub = p_ub->Predict(calc).FakeData(uboonePOT);

  SingleSampleExperiment expt_nd(p_nd, data_nd);
  SingleSampleExperiment expt_fd(p_fd, data_fd);
  SingleSampleExperiment expt_ub(p_ub, data_ub);

  MultiExperimentSBN multiExpt({&expt_nd, &expt_fd, &expt_ub}, {kSBND, kICARUS, kMicroBoone});
  MultiExperimentSBN fd_nd({&expt_nd, &expt_fd}, {kSBND, kICARUS});

  Surface surf_nom(&multiExpt, calc, kAxForTh, kAxDmSq);
  Surface surf_nd_fd(&fd_nd, calc, kAxForTh, kAxDmSq);
  Surface surf_nom_nd(&expt_nd, calc, kAxForTh, kAxDmSq);
  Surface surf_nom_fd(&expt_fd, calc, kAxForTh, kAxDmSq);
  Surface surf_nom_ub(&expt_ub, calc, kAxForTh, kAxDmSq);

  fout.mkdir("exclusion");
  fout.cd("exclusion");    

  surf_nom.SaveTo(gDirectory->mkdir("nom"));
  surf_nom_nd.SaveTo(gDirectory->mkdir("nom_nd"));
  surf_nom_fd.SaveTo(gDirectory->mkdir("nom_fd"));
  surf_nom_ub.SaveTo(gDirectory->mkdir("nom_ub"));
  surf_nd_fd.SaveTo(gDirectory->mkdir("nom_nd_fd"));

  std::map<std::string, std::vector<const ISyst*>> slists;
  std::vector<const ISyst*> systs = genie_systs;
  systs.insert(systs.end(), flux_systs.begin(), flux_systs.end());
  systs.push_back(kEnergyScaleMuon);
  systs.push_back(kEnergyScaleMuonND);
  systs.push_back(kEnergyScaleHadron);
  systs.push_back(kEnergyScaleHadronND);
  slists["all_systs"] = systs; 

  for(auto syst_pair: slists){
    std::string suffix = syst_pair.first;
    std::vector<const ISyst*> slist = syst_pair.second;    

    Surface surf_syst(&multiExpt, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_nd_fd(&fd_nd, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_nd(&expt_nd, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_fd(&expt_fd, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_ub(&expt_ub, calc, kAxForTh, kAxDmSq, {}, slist);

    surf_syst_nd.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
    surf_syst_fd.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
    surf_syst_ub.SaveTo(gDirectory->mkdir(("ub_"+suffix).c_str()));
    surf_syst.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
    surf_syst_nd_fd.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));
  } // end for s


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
//  calc2->SetL(kBaselineSBND);
//  const Spectrum data_nd2 = p_nd->Predict(calc2).FakeData(sbndPOT);
//  calc2->SetL(kBaselineIcarus);
//  const Spectrum data_fd2 = p_fd->Predict(calc2).FakeData(icarusPOT);
//  calc2->SetL(kBaselineMicroBoone);
//  const Spectrum data_ub2 = p_ub->Predict(calc2).FakeData(uboonePOT);
//
//  SingleSampleExperiment expt_nd2(p_nd, data_nd2);
//  SingleSampleExperiment expt_fd2(p_fd, data_fd2);
//  SingleSampleExperiment expt_ub2(p_ub, data_ub2);
//
//  MultiExperimentSBN multiExpt2({&expt_nd2, &expt_fd2, &expt_ub2}, {kSBND, kICARUS, kMicroBoone});
//  MultiExperimentSBN fd_nd2({&expt_nd2, &expt_fd2}, {kSBND, kICARUS});
//
//  //Surface surf_nom2(&multiExpt2, calc2, kAxForTh, kAxDmSq);
//  //Surface surf_nd_fd2(&fd_nd2, calc2, kAxForTh, kAxDmSq);
//  calc2->SetL(kBaselineSBND);
//  //Surface surf_nom_nd2(&expt_nd2, calc2, kAxForTh, kAxDmSq);
//  calc2->SetL(kBaselineIcarus);
//  //Surface surf_nom_fd2(&expt_fd2, calc2, kAxForTh, kAxDmSq);
//  calc2->SetL(kBaselineMicroBoone);
//  //Surface surf_nom_ub2(&expt_ub2, calc2, kAxForTh, kAxDmSq);
//    
//  //fout.cd("..");
//  //fout.mkdir("allowed");
//  //fout.cd("allowed");
//
//  //surf_nom2.SaveTo(gDirectory->mkdir("nom"));
//  //surf_nom_nd2.SaveTo(gDirectory->mkdir("nom_nd"));
//  //surf_nom_fd2.SaveTo(gDirectory->mkdir("nom_fd"));
//  //surf_nom_ub2.SaveTo(gDirectory->mkdir("nom_ub"));
//  //surf_nd_fd2.SaveTo(gDirectory->mkdir("nom_nd_fd"));
//
//  for(const std::vector<const ISyst*> slist: slists){
//    //Surface surf_syst2(&multiExpt2, calc2, kAxForTh, kAxDmSq, {}, slist);
//    //Surface surf_syst_nd_fd2(&fd_nd2, calc2, kAxForTh, kAxDmSq, {}, slist);
//    calc2->SetL(kBaselineSBND);
//    //Surface surf_syst_nd2(&expt_nd2, calc2, kAxForTh, kAxDmSq, {}, slist);
//    calc2->SetL(kBaselineIcarus);
//    //Surface surf_syst_fd2(&expt_fd2, calc2, kAxForTh, kAxDmSq, {}, slist);
//    calc2->SetL(kBaselineMicroBoone);
//    //Surface surf_syst_ub2(&expt_ub2, calc2, kAxForTh, kAxDmSq, {}, slist);
//
//    std::string suffix = "prop_systs";
//
//    //surf_syst_nd2.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
//    //surf_syst_fd2.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
//    //surf_syst_ub2.SaveTo(gDirectory->mkdir(("ub_"+suffix).c_str()));
//    //surf_syst2.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
//    //surf_syst_nd_fd2.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));
//
//  } // end for s
}
