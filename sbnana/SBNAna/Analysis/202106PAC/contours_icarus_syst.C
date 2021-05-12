#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/Prediction/PredictionInterp.h"
#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/Analysis/FitAxis.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Vars/FitVarsSterileApprox.h"

#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"

#include "sbnana/CAFAna/Analysis/Surface.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/SBNFluxSysts.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"
#include "TRandom3.h"

using namespace ana;

// NB these have locked in all the selection cuts
const std::string icarus_wildcard = "/pnfs/sbn/persistent/analysis/CAF/202106PAC/workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus/combined_reduced_flatcafs/*.flat.caf.root";

void contours_icarus_syst()
{
  Loaders loaders;
  loaders.SetLoaderPath(icarus_wildcard, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);

  const Var kTrueE = SIMPLEVAR(truth.E);

  const Var kRecoE([](const caf::SRSliceProxy* sr)
                   {
                     double E = sr->truth.E;
                     if(kIsNC(sr)) E *= sr->truth.inelasticityY;
                     return E * gRandom->Gaus(1, .1);
                   });

  const Var kCryostatWeight = Constant(2); // only using cryo0

  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSel = kNuMuCC_FullSelection;

  // TODO should make this an official function somewhere
  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                                   0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                   1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);


  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  const std::vector<const ISyst*> xsec_systs = GetSBNGenieWeightSysts();
  const std::vector<const ISyst*> flux_systs = GetSBNFluxHadronSysts(10);

  std::vector<const ISyst*> all_systs = xsec_systs;
  all_systs.insert(all_systs.end(), flux_systs.begin(), flux_systs.end());

  osc::NoOscillations noosc;

  NoExtrapPredictionGenerator gen(axEnergy, kNumuSpillSel, kNumuSel);
  PredictionInterp pred(all_systs, &noosc, gen, loaders);

  loaders.Go();

  // Fake data spectrum
  const Spectrum fake = pred.Predict(&noosc).FakeData(kPOTnominal);

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  const FitAxis kAxTh(&kFitSinSq2ThetaMuMu, 30, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 30, 1e-2, 1e2, true);
  
  SingleSampleExperiment expt(&pred, fake);

  Surface surf(&expt, calc, kAxTh, kAxDmSq, {}, all_systs);

  TH2* crit5sig = Gaussian5Sigma1D1Sided(surf);
  surf.DrawContour(crit5sig, kSolid, kBlack);

  gPad->Print("contour_icarus_syst.pdf");
}
