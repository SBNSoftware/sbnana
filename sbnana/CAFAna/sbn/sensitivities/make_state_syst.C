#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Analysis/Calcs.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Prediction/PredictionInterp.h"
#include "sbnana/CAFAna/Prediction/PredictionGenerator.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/SBNFluxSysts.h"
#include "sbnana/CAFAna/Systs/EnergySysts.h"
#include "sbnana/CAFAna/Systs/Systs.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "OscLib/IOscCalc.h"

#include "TFile.h"

#include <string>

using namespace ana;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void make_state_syst(const std::string anatype = numuStr)
{
  Loaders loaders_nd, loaders_fd, loaders_ub;

  if(anatype == numuStr) {
    const std::string dir = "/sbnd/data/users/jlarkin/cafs/sbnd/";
    const std::string fnameBeam_nd = dir + "sbnd_numu.flat.root";
    //const std::string fnameBeam_fd = dir + "output_SBNOsc_NumuSelection_Modern_Icarus.flat.root";
    //const std::string fnameBeam_ub = dir + "output_SBNOsc_NumuSelection_Modern_Uboone.flat.root";
    std::vector<std::string> fnameBeam_fd;
    //fnameBeam_fd.push_back("/pnfs/sbnd/scratch/users/jlarkin/cafs/icarus/reduced_icarus_1.flat.root");
    //fnameBeam_fd.push_back("reduced_icarus_0.flat.root");
    fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_1.flat.root");
    fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_2.flat.root");
    fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_3.flat.root");
    fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_4.flat.root");
    fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_5.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_6.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_7.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_8.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_9.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_10.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_11.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_12.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_13.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_14.flat.root");
    //fnameBeam_fd.push_back("/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_15.flat.root");

    std::string fnameBeam_fd_wildcard = "/pnfs/sbnd/persistent/users/jlarkin/cafs/icarus/icarus_numu_*.flat.root";

    loaders_nd.SetLoaderPath(fnameBeam_nd, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loaders_fd.SetLoaderPath(fnameBeam_fd_wildcard, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    //loaders_fd.SetLoaderFiles(fnameBeam_fd, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    //loaders_ub.SetLoaderPath(fnameBeam_ub, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  }

  else if (anatype == nueStr) {
    const std::string dir = "/sbnd/data/users/bzamoran/Modern_6thNov2019/";

    //BNB files contain nominal non-swap beam (so numubg, nuebg, NC)
    const std::string fnameBNB_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_Numu.flat.root";
    const std::string fnameBNB_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_Numu.flat.root";
    const std::string fnameBNB_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_Numu.flat.root";

    //Nue instrinsic only (to increase stats)
    const std::string fnameInt_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_Int.flat.root";
    const std::string fnameInt_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_Int.flat.root";
    const std::string fnameInt_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_Int.flat.root";

    //Swap files are for signal
    const std::string fnameSwap_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_Osc.flat.root";
    const std::string fnameSwap_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_Osc.flat.root";
    const std::string fnameSwap_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_Osc.flat.root";

    // //Dirt (background)
    // const std::string fnameDirt_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_CosDirt.flat.root";
    // const std::string fnameDirt_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_CosDirt.flat.root";
    // const std::string fnameDirt_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_CosDirt.flat.root";

    loaders_nd.SetLoaderPath( fnameBNB_nd,   Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders_nd.SetLoaderPath( fnameInt_nd,   Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders_nd.SetLoaderPath( fnameSwap_nd,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

    loaders_fd.SetLoaderPath( fnameBNB_fd,   Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders_fd.SetLoaderPath( fnameInt_fd,   Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders_fd.SetLoaderPath( fnameSwap_fd,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

    loaders_ub.SetLoaderPath( fnameBNB_ub,   Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders_ub.SetLoaderPath( fnameInt_ub,   Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders_ub.SetLoaderPath( fnameSwap_ub,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);
  }
  else {
    std::cout << "Unrecognised analysis - use numu or nue" << std::endl;
    return;
  }

  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kRecoE = SIMPLEVAR(freco.nuE);
  const Var kWeight = SIMPLEVAR(freco.wgt);
  const Var kWeightUB = Scaled(kWeight, 2.0*80.0/500.0);
  const Var kWeightFD = Scaled(kWeight, 2.0);
  //const Var kWeightFD = Constant(2.0);
  const Var kTrueWeight([](const caf::SRSliceProxy* slc)
         {
           if (abs(slc->truth.pdg) == 14 && slc->truth.iscc) return 0.8;
           return 1.0;
         });

  const Cut kFakeRecoMatched([](const caf::SRSliceProxy* slc)
         {
           return kHasMatchedNu(slc) && !std::isnan(slc->freco.nuE)
                                     && !std::isinf(slc->freco.nuE);
         });

  const Cut kNumuSelFD = kNuMuCC_FullSelection;
  const Cut kNumuSelND = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut;

  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                           0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                           1.25, 1.5, 2.0, 2.5, 3.0};
  const vector<double> fewer_binEdges = {0.3, 0.4, 0.45, 0.5,
                           0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                           1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  NoExtrapPredictionGenerator nom_gen(axEnergy, kNoSpillCut, kFakeRecoMatched, kWeight);
  //NoExtrapPredictionGenerator nom_gen(axEnergy, kNoSpillCut, kNumuSelND);
  //NoExtrapPredictionGenerator nom_genFD(axEnergy, kNoSpillCut, kNumuSelFD, kWeightFD);
  NoExtrapPredictionGenerator nom_genFD(axEnergy, kNoSpillCut, kFakeRecoMatched, kWeightFD);
  NoExtrapPredictionGenerator nom_genUB(axEnergy, kNoSpillCut, kFakeRecoMatched, kWeightUB);
  //NoExtrapPredictionGenerator nom_gen(axEnergy, kNoSpillCut, kHasMatchedNu);

  std::vector<const ISyst*> systs = GetSBNGenieWeightSysts();
  for(const ISyst* s: GetSBNFluxHadronSysts(30)) systs.push_back(s);
  for(const ISyst* s: GetEnergySysts()) systs.push_back(s);
  for(const ISyst* s: GetBigEnergySysts()) systs.push_back(s);
  //for(const ISyst* s: GetMiscSysts()) systs.push_back(s);
  systs.push_back(&GetMECSyst());
  systs.push_back(&GetPOTSyst());
  systs.push_back(&GetNormSyst());
  systs.push_back(&GetNormSystFD());
  systs.push_back(&GetNormSystND());

  osc::NoOscillations calc;

  PredictionInterp pred_nd(systs, &calc, nom_gen, loaders_nd);
  PredictionInterp pred_fd(systs, &calc, nom_genFD, loaders_fd);
  PredictionInterp pred_ub(systs, &calc, nom_genUB, loaders_fd);

//  loaders_nd.Go();
  loaders_fd.Go();
//  loaders_ub.Go();

  std::cout << "Creating file " << ("cafe_state_syst_"+anatype+".root").c_str() << std::endl;

  TFile fout(("cafe_state_syst_"+anatype+"_uboone.root").c_str(), "RECREATE");

  //pred_nd.SaveTo(fout.mkdir("pred_nd"));
  //pred_fd.SaveTo(fout.mkdir("pred_fd"));
  pred_ub.SaveTo(fout.mkdir("pred_ub"));
}



