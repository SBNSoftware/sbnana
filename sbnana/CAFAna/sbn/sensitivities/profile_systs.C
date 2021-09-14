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
#include "sbnana/CAFAna/Systs/SBNFluxSysts.h"
#include "sbnana/CAFAna/Systs/EnergySysts.h"
#include "sbnana/CAFAna/Systs/Systs.h"
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

void profile_systs(const std::string anatype = numuStr)
{
  const auto genie_systs = GetSBNGenieWeightSysts();
  const auto flux_systs = GetSBNFluxHadronSysts(30);
  const auto energy_systs = GetEnergySysts();
  GetMECSyst();
  GetPOTSyst();
  GetNormSyst();
  GetNormSystND();
  GetNormSystFD();
  
  const char* name_in;
  const char* name_out;
  if (anatype == numuStr) {
    name_in = "cafe_state_syst_numu.root";
    name_out = "profile_numu.root";
  }
  else if (anatype == nueStr) {
    name_in = "cafe_state_syst_nue.root";
    name_out = "profile_nue.root";
  }
  else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }

  TFile fin(name_in);
  TFile fout(name_out,"RECREATE");

  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  std::map<std::string, std::vector<const ISyst*>> slists;
  std::vector<const ISyst*> systs = genie_systs;
  systs.insert(systs.end(), flux_systs.begin(), flux_systs.end());
  systs.insert(systs.end(), energy_systs.begin(), energy_systs.end());
  systs.push_back(&GetMECSyst());
  systs.push_back(&GetNormSyst());
  systs.push_back(&GetNormSystND());
  systs.push_back(&GetNormSystFD());

  osc::NoOscillations noosc;
  TCanvas c;
  c.Print("energy_syst_profile.pdf[", "pdf");
  for(auto syst: systs) {
    if(syst->ShortName().find("IntraNuke") == std::string::npos) continue;
    std::vector<const ISyst *> loop_systs;
    std::copy_if(systs.begin(), systs.end(), std::back_inserter(loop_systs), 
                 [syst](auto s){ return s->ShortName() != syst->ShortName(); });
    auto title = syst->ShortName() + ";Shift (#sigma);#chi^{2}";
    auto title_fd = syst->ShortName() + " (fd);Shift (#sigma);#chi^{2}";
    auto data_nd = p_nd->Predict(&noosc).FakeData(sbndPOT);
    auto data_fd = p_fd->Predict(&noosc).FakeData(icarusPOT);
    SingleSampleExperiment expt_nd(p_nd, data_nd);
    SingleSampleExperiment expt_fd(p_fd, data_fd);
    MultiExperimentSBN fd_nd({&expt_nd, &expt_fd}, {kSBND, kICARUS});
    double x_min = -3.05;
    double x_max = 3.05;
    int nbins = 61;
    auto prof = Profile(&fd_nd, calc, syst, nbins, x_min, x_max, -1, {}, loop_systs);
    auto prof_fd = Profile(&expt_fd, calc, syst, nbins, x_min, x_max, -1, {}, loop_systs);
    auto slice = Profile(&fd_nd, calc, syst, nbins, x_min, x_max, -1, {}, {});
    auto slice_fd = Profile(&expt_fd, calc, syst, nbins, x_min, x_max, -1, {}, {});

    prof->SetTitle((syst->ShortName() + " - With All Other Systs").c_str());
    prof_fd->SetLineColor(kMagenta);
    prof_fd->SetMarkerColor(kMagenta);
    prof->SetMarkerStyle(kFullCircle);
    prof_fd->SetMarkerStyle(kFullSquare);
    prof->Draw("hist lp");
    prof_fd->Draw("hist lp same");
    TLegend * lgdis = new TLegend(0.35,0.60,0.65,0.89);
    lgdis->SetFillColor(0);
    lgdis->SetBorderSize(0);
    lgdis->AddEntry(prof_fd, "ICARUS Only fit");
    lgdis->AddEntry(prof, "SBND+ICARUS fit");
    lgdis->Draw("same");
    c.Print("energy_syst_profile.pdf", "pdf");
    slice->SetTitle((syst->ShortName() + " - No Other Systs").c_str());
    slice_fd->SetLineColor(kMagenta);
    slice_fd->SetMarkerColor(kMagenta);
    slice->SetMarkerStyle(kFullCircle);
    slice_fd->SetMarkerStyle(kFullSquare);
    slice->Draw("hist lp");
    slice_fd->Draw("hist lp same");
    lgdis->Draw("same");
    c.Print("energy_syst_profile.pdf", "pdf");
  }
  c.Print("energy_syst_profile.pdf]", "pdf");
}
