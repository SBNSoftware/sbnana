#pragma once

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/PrimaryUtils.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

namespace ana{

  /// \ref Neutrino pdg
  const Var kNuMITruePDG([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return 0;
    return PrimaryUtil::NeutrinoPDG_True(slc->truth);
  });
  /// \ref Target pdg
  const Var kNuMITrueTarget([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return 0;
    return PrimaryUtil::Target_True(slc->truth);
  });
  /// \ref GENIE interaction code (https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360)
  const Var kNuMITrueMode([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return PrimaryUtil::NeutrinoMode_True(slc->truth);
  });
  /// \ref IsCC (0:NC, 1:CC, -1:Not neutrino)
  const Var kNuMITrueIsCC([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    if (slc->truth.iscc) return 1;
    else return 0;
  });
  /// \ref Number of primary pi+
  const Var kNuMITrueNpip([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return PrimaryUtil::Npip_True(slc->truth);
  });
  /// \ref Number of primary pi-
  const Var kNuMITrueNpim([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return PrimaryUtil::Npim_True(slc->truth);
  });
  /// \ref Number of primary pi0
  const Var kNuMITrueNpi0([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return PrimaryUtil::Npi0_True(slc->truth);
  });
  /// \ref Neutrino energy
  const Var kNuMITrueNuE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return PrimaryUtil::NeutrinoE_True(slc->truth);
  });

  /// \ref 0: not signal, 1: signal
  const Var kNuMIIsSignal([](const caf::SRSliceProxy* slc) -> int {
    if( kNuMI_1muNp0piStudy_Signal_NoContainment_ProtonThreshold(slc) ) return 1;
    else return 0;
  });

  /// \ref 0: Muon candidate track exiting, 1: Muon candidate track contained (-1: no muon candidate)
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
