#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  // SpillCut

  // - CRTPMT matching

  extern const SpillCut spillcutHasValidFlash;
  extern const SpillCut spillcutHasInTimeFlash;
  extern const SpillCut spillcutCRTPMTCosmicByID;
  extern const SpillCut spillcutCRTPMTHasNegativeTOF;
  extern const SpillCut spillcutCRTPMTAllNegativeTOF;

  // Cut (slice)

  // - FV

  extern const Cut cutRecoVertexFV;
  extern const Cut cutFMScore;
  extern const Cut cutFMTime;

}
