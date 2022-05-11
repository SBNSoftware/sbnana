#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "cafanacore/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/Analysis/FitAxis.h"

#include "sbnana/CAFAna/Vars/FitVarsSterileApprox.h"

#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"

#include "sbnana/CAFAna/Analysis/Surface.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "TGraph.h"
#include "TPad.h"

using namespace ana;

// Since I'm in the SBND group the SAM dataset works fine for me
const std::string sbnd_wildcard = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_flat_caf_sbnd";

void contours_sbnd()
{
  Loaders loaders;
  loaders.SetLoaderPath(sbnd_wildcard, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);

  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kRecoE = kTrueE; // TODO

  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSel = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut;

  // TODO should make this an official function somewhere
  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                                   0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                   1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);


  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  PredictionNoExtrap pred(loaders, axEnergy, kNumuSpillSel, kNumuSel);

  loaders.Go();


  // Fake data spectrum
  osc::NoOscillations noosc;
  const Spectrum fake = pred.Predict(&noosc).FakeData(kPOTnominal);

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  const FitAxis kAxTh(&kFitSinSq2ThetaMuMu, 30, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 30, 1e-2, 1e2, true);

  SingleSampleExperiment expt(&pred, fake);

  Surface surf(&expt, calc, kAxTh, kAxDmSq);

  TH2* crit5sig = Gaussian5Sigma1D1Sided(surf);
  surf.DrawContour(crit5sig, kSolid, kBlack);

  gPad->Print("contour_sbnd.pdf");

  TH2* crit3sig = Gaussian3Sigma1D1Sided(surf);
  TH2* crit90pc = Gaussian90Percent1D1Sided(surf);
  TH2* crit95pc = Gaussian95Percent1D1Sided(surf);
  TH2* crit99pc = Gaussian99Percent1D1Sided(surf);

  TFile fout("sbnd_numu_disapp_exclusion_pac2021.root", "RECREATE");
  surf.GetGraphs(crit3sig)[0]->Write("nd_stat_3sig");
  surf.GetGraphs(crit5sig)[0]->Write("nd_stat_5sig");
  surf.GetGraphs(crit90pc)[0]->Write("nd_stat_90pct");
  surf.GetGraphs(crit95pc)[0]->Write("nd_stat_95pct");
  surf.GetGraphs(crit99pc)[0]->Write("nd_stat_99pct");

  expt.SaveTo(fout, "nd_expt");
}
