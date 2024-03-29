#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/Prediction/PredictionSBNExtrap.h"
#include "sbnana/CAFAna/Prediction/PredictionInterp.h"
#include "sbnana/CAFAna/Analysis/Calcs.h"
#include "sbnana/CAFAna/Analysis/FitAxis.h"

#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Vars/FitVarsSterileApprox.h"

#include "sbnana/CAFAna/Analysis/Surface.h"
#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"
#include "sbnana/CAFAna/Experiment/MultiExperimentSBN.h"
#include "sbnana/CAFAna/Experiment/CountingExperiment.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/BoosterFluxSysts.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

using namespace ana;

// State file
const char* basicFname = "cafe_state_extrap_syst.root";

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;

void nus_extrap(const char* stateFname = basicFname)
{
  if (TFile(stateFname).IsZombie()){
    std::cout << "Run make_state_extrap_syst.C first!" << std::endl;
    return;
  }

  std::vector<const ISyst*> systs = GetSBNWeightSysts();

  // Make sure all hadronization systs exist
  for(const ISyst* s: GetBoosterFluxHadronSysts(30));
  // But only use the first 10 of them
  for(const ISyst* s: GetBoosterFluxHadronSysts(10)) systs.push_back(s);

  std::cout << "Loading state from " << stateFname << std::endl; 
  TFile fin(stateFname);
  IPrediction* pred = ana::LoadFrom<IPrediction>(fin.GetDirectory("pred_syst")).release();
  IPrediction* predND = ana::LoadFrom<IPrediction>(fin.GetDirectory("predND_syst")).release();

  // Calculator
  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();
  // TODO - why is this necessary?
  ((PredictionInterp*)pred)->SetOscSeed(calc);

  // To make a fit we need to have a "data" spectrum to compare to our MC
  // Prediction object
  calc->SetL(kBaselineIcarus);
  const Spectrum data = pred->Predict(calc).FakeData(icarusPOT);
  SingleSampleExperiment expt(pred, data);

  calc->SetL(kBaselineSBND);
  const Spectrum dataND = predND->Predict(calc).FakeData(sbndPOT);
  CountingExperiment exptCount(predND, dataND);

  MultiExperimentSBN multiExpt({&expt, &exptCount}, {kICARUS, kSBND});

  //Define fit axes
  const FitAxis kAxSinSq2ThetaMuMu(&kFitSinSq2ThetaMuMu, 40, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 40, 2e-2, 1e2, true);

  const FitAxis kAxSinSq2ThetaMuMuCoarse(&kFitSinSq2ThetaMuMu, 10, 1e-3, 1, true);
  const FitAxis kAxDmSqCoarse(&kFitDmSqSterile, 10, 2e-2, 1e2, true);

  TCanvas* c1 = new TCanvas("c1");
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.15);

  Surface surf(&expt, calc,
               kAxSinSq2ThetaMuMuCoarse,
               kAxDmSqCoarse,
               {}, systs);

  Surface surfStat(&expt, calc,
                   kAxSinSq2ThetaMuMu,
                   kAxDmSq);

  calc->SetL(kBaselineSBND);
  Surface surfCount(&exptCount, calc,
                    kAxSinSq2ThetaMuMuCoarse,
                    kAxDmSqCoarse,
                    {}, systs);

  Surface surfCountStat(&exptCount, calc,
                        kAxSinSq2ThetaMuMu,
                        kAxDmSq);

  Surface surfJoint(&multiExpt, calc,
                    kAxSinSq2ThetaMuMuCoarse,
                    kAxDmSqCoarse,
                    {}, systs);

  Surface surfJointStat(&multiExpt, calc,
                        kAxSinSq2ThetaMuMu,
                        kAxDmSq);

  TH2* crit3sig = Gaussian3Sigma1D1Sided(surfStat);
  TH2* crit3sigCoarse = Gaussian3Sigma1D1Sided(surf);

  surfJointStat.DrawContour(crit3sig, 7, kRed);
  surfJoint.DrawContour(crit3sigCoarse, kSolid, kRed);

  surfStat.DrawContour(crit3sig, 7, kGreen+2);
  surf.DrawContour(crit3sigCoarse, kSolid, kGreen+2);

  surfCountStat.DrawContour(crit3sig, 7, kBlue);
  surfCount.DrawContour(crit3sigCoarse, kSolid, kBlue);

  TH1F* hDummy1 = new TH1F("","",1,0,1);
  TH1F* hDummy2 = new TH1F("","",1,0,1);
  TH1F* hDummy3 = new TH1F("","",1,0,1);

  hDummy1->SetLineColor(kRed);
  hDummy2->SetLineColor(kBlue);
  hDummy3->SetLineColor(kGreen+2);

  TLegend* leg = new TLegend(0.16, 0.16, 0.6, 0.4);
  leg->SetFillStyle(0);
  leg->SetHeader("3#sigma C.L.","C");
  leg->AddEntry(hDummy2, "SBND counting only", "l");
  leg->AddEntry(hDummy3, "Icarus/SBND ratio only", "l");
  leg->AddEntry(hDummy1, "Combined fit", "l");
  leg->Draw();

  gPad->Print("nus_extrap.pdf");
}
