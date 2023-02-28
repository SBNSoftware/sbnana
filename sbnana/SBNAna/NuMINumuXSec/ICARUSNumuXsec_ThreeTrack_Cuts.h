#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_ThreeTrack_Variables.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace ThreeTrack{

  extern const Cut Pass3T1PCut;

  extern const Cut HasThreePrimaryTracks;
  extern const Cut HasMuonTrack;
  extern const Cut HasProtonTrack;
  extern const Cut MuonContained;
  extern const Cut MuonExiting;

} // end namespace ThreeTrack

} // end namespace ICARUSNumuXsec
