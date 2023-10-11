#include "sbnana/SBNAna/Cuts/NuMIXSecRecoEventVars.h"

namespace ana{

  // Neutrino pdg
  const Var kNuMITruePDG([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return 0;
    return kTruth_NeutrinoPDG(&slc->truth);
  });
  // Target pdg
  const Var kNuMITrueTarget([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return 0;
    return kTruth_Target(&slc->truth);
  });
  // GENIE interaction code (https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360)
  const Var kNuMITrueMode([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_NeutrinoMode(&slc->truth);
  });
  // IsCC (0:NC, 1:CC, -1:Not neutrino)
  const Var kNuMITrueIsCC([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    if (slc->truth.iscc) return 1;
    else return 0;
  });
  // Number of primary proton
  const Var kNuMITrueNProton([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_NProton_Primary(&slc->truth);
  });
  // Number of primary neutron
  const Var kNuMITrueNNeutron([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_NNeutron_Primary(&slc->truth);
  });
  // Number of primary pi+
  const Var kNuMITrueNpip([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npip_Primary(&slc->truth);
  });
  // Number of ANY pi+ (not just primary)
  const Var kNuMITrueNpip_All([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npip_All(&slc->truth);
  });
  // Number of primary pi-
  const Var kNuMITrueNpim([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npim_Primary(&slc->truth);
  });
  // Number of ANY pi- (not just primary)
  const Var kNuMITrueNpim_All([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npim_All(&slc->truth);
  });
  // Number of primary pi0
  const Var kNuMITrueNpi0([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npi0_Primary(&slc->truth);
  });
  // Number of ANY pi0 (not just primary)
  const Var kNuMITrueNpi0_All([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npi0_All(&slc->truth);
  });
  // Neutrino energy
  const Var kNuMITrueNuE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_NeutrinoE(&slc->truth);
  });
  // Q2
  const Var kNuMITrueQ2([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_Q2(&slc->truth);
  });
  // q0; energy transfer
  const Var kNuMITrueq0([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_q0(&slc->truth);
  });
  // q3; momentum transfer
  const Var kNuMITrueq3([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_q3(&slc->truth);
  });
  // w; hadronic mass
  const Var kNuMITruew([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_w(&slc->truth);
  });
  // 0: RHC, 1: FHC
  const Var kNuMIIsFHC([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -5; //TODO Define better dummy value
    return kTruth_IsFHC(&slc->truth);
  });

  // True muon kinetic energy
  const Var kNuMITrueMuonKE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_MuonKE(&slc->truth);
  });
  // True muon cosine angle w.r.t. neutrino
  const Var kNuMITrueMuonNuCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_MuonNuCosineTheta(&slc->truth);
  });
  // True muon contain?: 1: contained, 0: not contained (-1: muon not found)
  const Var kNuMITrueMuonContained([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1; //TODO Define better dummy value
    return kTruth_MuonContained(&slc->truth);
  });
  // True proton kinetic energy
  const Var kNuMITrueProtonKE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_ProtonKE(&slc->truth);
  });
  // True proton cosine angle w.r.t. neutrino
  const Var kNuMITrueProtonNuCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_ProtonNuCosineTheta(&slc->truth);
  });

  // True Charged pion kinetic energy
  const Var kNuMITrueChargedPionKE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_ChargedPionKE(&slc->truth);
  });

  // 0: Muon candidate track exiting, 1: Muon candidate track contained (-1: no muon candidate)
  const Var kNuMIRecoMuonContained([](const caf::SRSliceProxy* slc) -> int {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });

}