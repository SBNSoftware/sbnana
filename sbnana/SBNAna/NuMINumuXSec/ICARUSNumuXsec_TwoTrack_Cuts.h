#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Variables.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  extern const Cut HasTwoPrimaryTracks;
  extern const Cut HasMuonTrack;
  extern const Cut HasProtonTrack;
  extern const Cut MuonTrackContained;
  extern const Cut MuonTrackExiting;

  extern const Cut MuonTrackTruthContainedNuMuon;
  extern const Cut MuonTrackTruthExitingNuMuon;
  extern const Cut MuonTrackTruthCosmicMuon;
  extern const Cut MuonTrackTruthStoppingProton;
  extern const Cut MuonTrackTruthInelProton;
  extern const Cut MuonTrackTruthOtherProton;
  extern const Cut MuonTrackTruthStoppingChargedPion;
  extern const Cut MuonTrackTruthExitingChargedPion;
  extern const Cut MuonTrackTruthOther;

  extern const Cut ProtonPCut;

  // pion tagging
  extern const Cut NoPionCand;

  namespace Aux{
    extern const Cut HasRelaxedMuonTrack;
    extern const Cut HasRelaxedProtonTrack;
    extern const Cut RelaxedMuonTrackContained;
    extern const Cut RelaxedProtonTrackContained;

    extern const Cut RelaxedMuonTrackTruthContainedNuMuon;
    extern const Cut RelaxedMuonTrackTruthExitingNuMuon;
    extern const Cut RelaxedMuonTrackTruthCosmicMuon;
    extern const Cut RelaxedMuonTrackTruthStoppingProton;
    extern const Cut RelaxedMuonTrackTruthInelProton;
    extern const Cut RelaxedMuonTrackTruthOtherProton;
    extern const Cut RelaxedMuonTrackTruthStoppingChargedPion;
    extern const Cut RelaxedMuonTrackTruthExitingChargedPion;
    extern const Cut RelaxedMuonTrackTruthOther;

    extern const Cut RelaxedProtonTrackTruthContainedNuMuon;
    extern const Cut RelaxedProtonTrackTruthExitingNuMuon;
    extern const Cut RelaxedProtonTrackTruthCosmicMuon;
    extern const Cut RelaxedProtonTrackTruthStoppingProton;
    extern const Cut RelaxedProtonTrackTruthInelProton;
    extern const Cut RelaxedProtonTrackTruthOtherProton;
    extern const Cut RelaxedProtonTrackTruthStoppingChargedPion;
    extern const Cut RelaxedProtonTrackTruthExitingChargedPion;
    extern const Cut RelaxedProtonTrackTruthOther;

    extern const Cut RelaxedProtonTrackIsochronous;
    extern const Cut RelaxedProtonTrackHasChi2MuonExcess;
    extern const Cut RelaxedProtonTrackHasChi2MuonDeficit;
    extern const Cut RelaxedProtonTrackShort;

  }

} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
