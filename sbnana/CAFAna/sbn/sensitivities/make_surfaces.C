#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Vars/FitVarsSterileApprox.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/CountingExperiment.h"
#include "CAFAna/Analysis/ExpInfo.h"
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Systs/SBNWeightSysts.h"
#include "CAFAna/Systs/UniverseOracle.h"
using namespace ana;

#include "OscLib/IOscCalculator.h"

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

void make_surfaces()
{
  const std::vector<const ISyst*>& systs = GetSBNWeightSysts();

  auto systs_flux = GetSBNFluxWeightSysts();
  auto systs_genie = GetSBNGenieWeightSysts();
  
  std::vector<const ISyst*> systs_to_process;

  std::vector<std::string> syst_names{"expskin_FluxUnisim","horncurrent_FluxUnisim","kminus_PrimaryHadronNormalization","kplus_PrimaryHadronFeynmanScaling","kzero_PrimaryHadronSanfordWang","nucleoninexsec_FluxUnisim","nucleonqexsec_FluxUnisim","nucleontotxsec_FluxUnisim","piminus_PrimaryHadronSWCentralSplineVariation","pioninexsec_FluxUnisim","pionqexsec_FluxUnisim","piontotxsec_FluxUnisim","piplus_PrimaryHadronSWCentralSplineVariation","genie_ccresAxial_Genie","genie_ncresAxial_Genie","genie_qema_Genie","genie_NC_Genie","genie_NonResRvbarp1pi_Genie","genie_NonResRvbarp2pi_Genie","genie_NonResRvp1pi_Genie","genie_NonResRvp2pi_Genie","genie_NonResRvbarp1piAlt_Genie","genie_NonResRvbarp2piAlt_Genie","genie_NonResRvp1piAlt_Genie","genie_NonResRvp2piAlt_Genie"};

  for (auto s : systs) {
    for (auto n : syst_names) if (n == s->ShortName()) systs_to_process.push_back(s);
  }

  TFile fin("cafe_state_syst_numu.root");
  TFile fout("surfaces.root","RECREATE");

  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
  PredictionInterp* p_ub = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_ub")).release();

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  //Define fit axes
  const FitAxis kAxSinSq2ThetaMuMu(&kFitSinSq2ThetaMuMu, 40, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 40, 1e-2, 1e2, true);

  // We'll call zero nominal
  calc->SetL(kBaselineSBND);
  const Spectrum data_nd = p_nd->Predict(calc).FakeData(sbndPOT);
  calc->SetL(kBaselineIcarus);
  const Spectrum data_fd = p_fd->Predict(calc).FakeData(icarusPOT);
  calc->SetL(kBaselineMicroBoone);
  const Spectrum data_ub = p_ub->Predict(calc).FakeData(uboonePOT);

  SingleSampleExperiment expt_nd(p_nd, data_nd);
  SingleSampleExperiment expt_fd(p_fd, data_fd);
  SingleSampleExperiment expt_ub(p_ub, data_ub);

  MultiExperimentSBN multiExpt({&expt_nd, &expt_fd, &expt_ub}, {kSBND, kICARUS, kMicroBoone});
  MultiExperimentSBN fd_nd({&expt_nd, &expt_fd}, {kSBND, kICARUS});

  Surface surf_nom(&multiExpt, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
  Surface surf_nd_fd(&fd_nd, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
  calc->SetL(kBaselineSBND);
  Surface surf_nom_nd(&expt_nd, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
  calc->SetL(kBaselineIcarus);
  Surface surf_nom_fd(&expt_fd, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
  calc->SetL(kBaselineMicroBoone);
  Surface surf_nom_ub(&expt_ub, calc, kAxSinSq2ThetaMuMu, kAxDmSq);

  fout.mkdir("exclusion");
  fout.cd("exclusion");    

  surf_nom.SaveTo(gDirectory->mkdir("nom"));
  surf_nom_nd.SaveTo(gDirectory->mkdir("nom_nd"));
  surf_nom_fd.SaveTo(gDirectory->mkdir("nom_fd"));
  surf_nom_ub.SaveTo(gDirectory->mkdir("nom_ub"));
  surf_nd_fd.SaveTo(gDirectory->mkdir("nom_nd_fd"));

  std::vector<std::vector<const ISyst*>> slists;
  //slists.emplace_back(1, systs[0]);
  slists.push_back(systs_to_process);

  for(const std::vector<const ISyst*> slist: slists){
    Surface surf_syst(&multiExpt, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);
    Surface surf_syst_nd_fd(&fd_nd, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);
    calc->SetL(kBaselineSBND);
    Surface surf_syst_nd(&expt_nd, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);
    calc->SetL(kBaselineIcarus);
    Surface surf_syst_fd(&expt_fd, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);
    calc->SetL(kBaselineMicroBoone);
    Surface surf_syst_ub(&expt_ub, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);

    std::string suffix = "prop_systs";

    surf_syst_nd.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
    surf_syst_fd.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
    surf_syst_ub.SaveTo(gDirectory->mkdir(("ub_"+suffix).c_str()));
    surf_syst.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
    surf_syst_nd_fd.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));

  } // end for s


  // Allowed Region

  OscCalcSterileApproxAdjustable* calc2 = DefaultSterileApproxCalc();
  calc2->calc.SetSinSq2ThetaMuMu(4*0.135*0.135*(1-0.135*0.135));
  calc2->calc.SetDmsq(1.32);

  //Define fit axes
  const FitAxis kAxSinSq2ThetaMuMu2(&kFitSinSq2ThetaMuMu, 40, 1e-3, 1, true);
  const FitAxis kAxDmSq2(&kFitDmSqSterile, 40, 1e-1, 1e1, true);

  // We'll call zero nominal
  calc2->SetL(kBaselineSBND);
  const Spectrum data_nd2 = p_nd->Predict(calc2).FakeData(sbndPOT);
  calc2->SetL(kBaselineIcarus);
  const Spectrum data_fd2 = p_fd->Predict(calc2).FakeData(icarusPOT);
  calc2->SetL(kBaselineMicroBoone);
  const Spectrum data_ub2 = p_ub->Predict(calc2).FakeData(uboonePOT);

  SingleSampleExperiment expt_nd2(p_nd, data_nd2);
  SingleSampleExperiment expt_fd2(p_fd, data_fd2);
  SingleSampleExperiment expt_ub2(p_ub, data_ub2);

  MultiExperimentSBN multiExpt2({&expt_nd2, &expt_fd2, &expt_ub2}, {kSBND, kICARUS, kMicroBoone});
  MultiExperimentSBN fd_nd2({&expt_nd2, &expt_fd2}, {kSBND, kICARUS});

  Surface surf_nom2(&multiExpt2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2);
  Surface surf_nd_fd2(&fd_nd2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2);
  calc2->SetL(kBaselineSBND);
  Surface surf_nom_nd2(&expt_nd2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2);
  calc2->SetL(kBaselineIcarus);
  Surface surf_nom_fd2(&expt_fd2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2);
  calc2->SetL(kBaselineMicroBoone);
  Surface surf_nom_ub2(&expt_ub2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2);
    
  fout.cd("..");
  fout.mkdir("allowed");
  fout.cd("allowed");

  surf_nom2.SaveTo(gDirectory->mkdir("nom"));
  surf_nom_nd2.SaveTo(gDirectory->mkdir("nom_nd"));
  surf_nom_fd2.SaveTo(gDirectory->mkdir("nom_fd"));
  surf_nom_ub2.SaveTo(gDirectory->mkdir("nom_ub"));
  surf_nd_fd2.SaveTo(gDirectory->mkdir("nom_nd_fd"));

  for(const std::vector<const ISyst*> slist: slists){
    Surface surf_syst2(&multiExpt2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2, {}, slist);
    Surface surf_syst_nd_fd2(&fd_nd2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2, {}, slist);
    calc2->SetL(kBaselineSBND);
    Surface surf_syst_nd2(&expt_nd2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2, {}, slist);
    calc2->SetL(kBaselineIcarus);
    Surface surf_syst_fd2(&expt_fd2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2, {}, slist);
    calc2->SetL(kBaselineMicroBoone);
    Surface surf_syst_ub2(&expt_ub2, calc2, kAxSinSq2ThetaMuMu2, kAxDmSq2, {}, slist);

    std::string suffix = "prop_systs";

    surf_syst_nd2.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
    surf_syst_fd2.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
    surf_syst_ub2.SaveTo(gDirectory->mkdir(("ub_"+suffix).c_str()));
    surf_syst2.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
    surf_syst_nd_fd2.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));

  } // end for s
}