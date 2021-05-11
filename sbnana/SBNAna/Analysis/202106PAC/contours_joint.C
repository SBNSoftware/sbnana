#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/Experiment/MultiExperiment.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/Analysis/FitAxis.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Vars/FitVarsSterileApprox.h"

#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"

#include "sbnana/CAFAna/Analysis/Surface.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"
#include "TRandom3.h"

using namespace ana;

// This should be
// workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus
// but that is only available from SAM to icarus users, for now.
//const std::string icarus_wildcard = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/*/*/*.root";

// That dataset takes absolutely forever, because it isn't flattened. Use this
// file for now.
const std::string icarus_wildcard = "/icarus/data/users/mueller/NuMuSelection/v09_17_00/nucosmics/icarus_nucosmics.flat.root";


// Since I'm in the SBND group the SAM dataset works fine for me
const std::string sbnd_wildcard = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_flat_caf_sbnd";

void contours_joint()
{
  Loaders loadersND, loadersFD;
  loadersND.SetLoaderPath(sbnd_wildcard, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loadersFD.SetLoaderPath(icarus_wildcard, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);

  const Var kTrueE = SIMPLEVAR(truth.E);

  const Var kRecoE([](const caf::SRSliceProxy* sr)
                   {
                     double E = sr->truth.E;
                     if(kIsNC(sr)) E *= sr->truth.inelasticityY;
                     return E * gRandom->Gaus(1, .1);
                   });
  //  const Var kRecoE = kTrueE; // TODO

  const Var kCryostatWeight = Constant(2); // only using cryo0

  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSelND = kSlcIsRecoNu && kSlcNuScoreCut && kFiducialVolumeND && kSlcFlashMatchCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut;
  const Cut kNumuSelFD = kNuMuCC_FullSelection;

  // TODO should make this an official function somewhere
  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                                   0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                   1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);


  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);
  const HistAxis axTrueEnergy("True energy (GeV)", binsEnergy, kTrueE);

  PredictionNoExtrap predND(loadersND, axEnergy, kNumuSpillSel, kNumuSelND);
  PredictionNoExtrap predNDTrue(loadersND, axTrueEnergy, kNumuSpillSel, kNumuSelND);

  PredictionNoExtrap predFD(loadersFD, axEnergy, kNumuSpillSel, kNumuSelFD, kNoShift, kCryostatWeight);
  PredictionNoExtrap predFDTrue(loadersFD, axTrueEnergy, kNumuSpillSel, kNumuSelFD, kNoShift, kCryostatWeight);

  loadersND.Go();
  loadersFD.Go();

  // Fake data spectrum
  osc::NoOscillations noosc;
  const Spectrum fakeND = predND.Predict(&noosc).FakeData(kPOTnominal);
  const Spectrum fakeNDTrue = predNDTrue.Predict(&noosc).FakeData(kPOTnominal);
  const Spectrum fakeFD = predFD.Predict(&noosc).FakeData(kPOTnominal);
  const Spectrum fakeFDTrue = predFDTrue.Predict(&noosc).FakeData(kPOTnominal);

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  const FitAxis kAxTh(&kFitSinSq2ThetaMuMu, 30, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 30, 1e-2, 1e2, true);
  
  SingleSampleExperiment exptND(&predND, fakeND);
  SingleSampleExperiment exptNDTrue(&predNDTrue, fakeNDTrue);
  SingleSampleExperiment exptFD(&predFD, fakeFD);
  SingleSampleExperiment exptFDTrue(&predFDTrue, fakeFDTrue);

  MultiExperiment multi({&exptND, &exptFD});
  MultiExperiment multiTrue({&exptNDTrue, &exptFDTrue});

  Surface surfND(&exptND, calc, kAxTh, kAxDmSq);
  Surface surfNDTrue(&exptNDTrue, calc, kAxTh, kAxDmSq);
  Surface surfFD(&exptFD, calc, kAxTh, kAxDmSq);
  Surface surfFDTrue(&exptFDTrue, calc, kAxTh, kAxDmSq);
  Surface surfMulti(&multi, calc, kAxTh, kAxDmSq);
  Surface surfMultiTrue(&multiTrue, calc, kAxTh, kAxDmSq);

  TH2* crit5sig = Gaussian5Sigma1D1Sided(surfND);
  surfND.DrawContour(crit5sig, kSolid, kBlack);
  surfNDTrue.DrawContour(crit5sig, kSolid, kRed);
  surfFD.DrawContour(crit5sig, kSolid, kBlack);
  surfFDTrue.DrawContour(crit5sig, kSolid, kRed);
  surfMulti.DrawContour(crit5sig, kSolid, kBlack);
  surfMultiTrue.DrawContour(crit5sig, kSolid, kRed);

  gPad->Print("contour_joint.pdf");
  gPad->Print("contour_joint.png");
}
