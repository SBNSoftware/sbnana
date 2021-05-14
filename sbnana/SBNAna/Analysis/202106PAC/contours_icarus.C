#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
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

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"

#include "TPad.h"

using namespace ana;

// This should be
// workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus
// but that is only available from SAM to icarus users, for now.
//const std::string icarus_wildcard = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/*/*/*.root";

// That dataset takes absolutely forever, because it isn't flattened. Use this
// file for now.
const std::string icarus_wildcard = "/icarus/data/users/mueller/NuMuSelection/v09_17_00/nucosmics/icarus_nucosmics.flat.root";

void contours_icarus()
{
  Loaders loaders;
  loaders.SetLoaderPath(icarus_wildcard, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);

  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kRecoE = kTrueE; // TODO

  const Var kCryostatWeight = Constant(2); // only using cryo0

  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSel = kNuMuCC_FullSelection;

  // TODO should make this an official function somewhere
  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                                   0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                   1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);


  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  PredictionNoExtrap pred(loaders, axEnergy, kNumuSpillSel, kNumuSel, kNoShift, kCryostatWeight);

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

  gPad->Print("contour_icarus.pdf");
}
