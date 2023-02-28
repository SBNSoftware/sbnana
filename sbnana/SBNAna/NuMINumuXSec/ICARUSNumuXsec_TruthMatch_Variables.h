#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TruthMatch{

  // For a given true muon (truth_index), find a reco track whose best-matched is this particle
  extern const Var varTruthMuonIndex;
  extern const Var varTruthMuonMatchedTrackIndex;
  extern const Var varTruthMuonMatchedTrackContainedness;
  extern const Var varTruthMuonMatchedTrackChi2Proton;
  extern const Var varTruthMuonMatchedTrackChi2Muon;
  extern const Var varTruthMuonMatchedTrackChi2Pion;

  // For a given true proton (truth_index), find a reco track whose best-matched is this particle
  extern const Var varTruthProtonIndex;
  extern const Var varTruthProtonMatchedTrackIndex;
  extern const Var varTruthProtonMatchedTrackContainedness;
  extern const Var varTruthProtonMatchedTrackChi2Proton;
  extern const Var varTruthProtonMatchedTrackChi2Muon;
  extern const Var varTruthProtonMatchedTrackChi2Pion;

  // For a given true charged pion (truth_index), find a reco track whose best-matched is this particle
  extern const Var varTruthChargedPionIndex;
  extern const Var varTruthChargedPionMatchedTrackIndex;
  extern const Var varTruthChargedPionMatchedTrackContainedness;
  extern const Var varTruthChargedPionMatchedTrackEndProcess;


} // end namespace TruthMatch

} // end namespace ICARUSNumuXsec
