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
  extern const Var TruthMuonIndex;
  extern const Var TruthMuonLength;
  extern const Var TruthMuonKE;
  extern const Var TruthMuonMatchedTrackIndex;
  extern const Var TruthMuonMatchedTrackContainedness;
  extern const Var TruthMuonMatchedTrackChi2Proton;
  extern const Var TruthMuonMatchedTrackChi2Muon;
  extern const Var TruthMuonMatchedTrackChi2Pion;

  // For a given true proton (truth_index), find a reco track whose best-matched is this particle
  extern const Var TruthProtonIndex;
  extern const Var TruthProtonLength;
  extern const Var TruthProtonKE;
  extern const Var TruthProtonMatchedTrackIndex;
  extern const Var TruthProtonMatchedTrackContainedness;
  extern const Var TruthProtonMatchedTrackChi2Proton;
  extern const Var TruthProtonMatchedTrackChi2Muon;
  extern const Var TruthProtonMatchedTrackChi2Pion;

  // For a given true charged pion (truth_index), find a reco track whose best-matched is this particle
  extern const Var TruthChargedPionIndex;
  extern const Var TruthChargedPionLength;
  extern const Var TruthChargedPionKE;
  extern const Var TruthChargedPionMatchedTrackIndex;
  extern const Var TruthChargedPionMatchedTrackContainedness;
  extern const Var TruthChargedPionMatchedTrackEndProcess;


} // end namespace TruthMatch

} // end namespace ICARUSNumuXsec
