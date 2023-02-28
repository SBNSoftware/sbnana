#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Variables.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  extern const Cut HasTwoPrimaryTracks;
  extern const Cut HasMuonTrack;
  extern const Cut HasProtonTrack;
  extern const Cut MuonContained;
  extern const Cut MuonExiting;

  extern const Cut ProtonPCut;

  // pion tagging
  extern const Cut NoPionCand;

  namespace Aux{
    extern const Cut HasRelaxedMuonTrack;
    extern const Cut HasRelaxedProtonTrack;

    extern const Cut RelaxedMuonTrackTruthContainedNuMuon;
    extern const Cut RelaxedMuonTrackTruthExitingNuMuon;
    extern const Cut RelaxedMuonTrackTruthCosmicMuon;
    extern const Cut RelaxedMuonTrackTruthStoppingProton;
    extern const Cut RelaxedMuonTrackTruthOtherProton;
    extern const Cut RelaxedMuonTrackTruthOther;

    extern const Cut RelaxedProtonTrackTruthContainedNuMuon;
    extern const Cut RelaxedProtonTrackTruthExitingNuMuon;
    extern const Cut RelaxedProtonTrackTruthCosmicMuon;
    extern const Cut RelaxedProtonTrackTruthStoppingProton;
    extern const Cut RelaxedProtonTrackTruthOtherProton;
    extern const Cut RelaxedProtonTrackTruthOther;

  }

} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
