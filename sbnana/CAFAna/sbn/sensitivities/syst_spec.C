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
using namespace ana;

#include "OscLib/IOscCalc.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"
#include "TRandom.h"

#include <vector>

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void syst_spec(const std::string anatype = numuStr)
{
  std::map<std::string, std::vector<const ISyst*>> slists;
  slists["xsec"] = GetSBNGenieWeightSysts();
  slists["xsec"].push_back(&GetMECSyst());
  slists["flux"] = GetSBNFluxHadronSysts(30);
  slists["detector"] = GetEnergySysts();
  slists["detector"].push_back(&GetPOTSyst());
  slists["detector"].push_back(&GetNormSystND());
  slists["detector"].push_back(&GetNormSystFD());
  slists["all"] = slists["detector"];
  slists["all"].insert(slists["all"].end(), slists["xsec"].begin(), slists["xsec"].end());
  slists["all"].insert(slists["all"].end(), slists["flux"].begin(), slists["flux"].end());
  GetNormSyst();

  const char* name_in;
  const char* name_out;
  if (anatype == numuStr) {
    name_in = "cafe_state_syst_numu.root";
    name_out = "output_syst_spec_numu.root";
  } else if (anatype == nueStr) {
    name_in = "cafe_state_syst_nue.root";
    name_out = "output_syst_spec_nue.root";
  } else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }

  TFile fin(name_in);
  TFile fout(name_out,"RECREATE");

  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
 
  osc::NoOscillations unosc;

  int i = 0;

  p_nd->Predict(&unosc).ToTH1(sbndPOT)->Write("spec_nd_nominal");
  p_fd->Predict(&unosc).ToTH1(icarusPOT)->Write("spec_fd_nominal");
  
  for(auto [name, systs]: slists) {
    i++;
    std::vector<TH1*> hists;
    for (int j = 0; j < 1000; ++j){
      hists.push_back(p_nd->PredictSyst(&unosc, SystShifts::RandomThrow(systs)).ToTH1(sbndPOT));
    }
    double xbins[] = {0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.25, 1.5, 2., 2.5, 3.};
    TH1D *shifted = new TH1D(("h"+to_string(i)).c_str(), "hist",19,xbins);
    TH1D *lower = new TH1D(("h2"+to_string(i)).c_str(), "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      shifted->SetBinContent(k, bincont[840]);
      lower->SetBinContent(k, bincont[160]);
    }    
    
    shifted->Write(("spec_nd_"+name+"_+1").c_str(), TObject::kOverwrite);
    lower->Write(("spec_nd_"+name+"_-1").c_str(), TObject::kOverwrite);

    i++;

    std::vector<TH1*> hists2;
    for (int j = 0; j < 1000; ++j){
      hists2.push_back(p_fd->PredictSyst(&unosc, SystShifts::RandomThrow(systs)).ToTH1(icarusPOT));
    }
    TH1D *shifted2 = new TH1D(("h"+to_string(i)).c_str(), "hist",19,xbins);
    TH1D *lower2 = new TH1D(("h2"+to_string(i)).c_str(), "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists2) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      shifted2->SetBinContent(k, bincont[840]);
      lower2->SetBinContent(k, bincont[160]);
    }    
    
    shifted2->Write(("spec_fd_"+name+"_+1").c_str(), TObject::kOverwrite);
    lower2->Write(("spec_fd_"+name+"_-1").c_str(), TObject::kOverwrite);
  }
}
