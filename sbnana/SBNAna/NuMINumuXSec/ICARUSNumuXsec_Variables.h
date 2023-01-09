#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  // SpillVar
  extern const SpillVar spillvarCountSpill;
  // - Test
  extern const SpillVar spillvarTest;
  // - PMT-CRT matching
  extern const SpillMultiVar spillvarOpFlashTime;
  extern const SpillMultiVar spillvarValidOpFlashTime;
  extern const SpillMultiVar spillvarInTimeOpFlashTime; // Valid AND InTime
  extern const SpillMultiVar spillvarCRTPMTTime;
  extern const SpillMultiVar spillvarCRTPMTMatchingID;
  extern const SpillVar spillvarCRTPMTMatchingEventID;

  // Var
  // - Flash matching
  extern const Var varFMScore;
  extern const Var varFMTime;
  // - Longest track
  //   - index
  extern const Var varLongestTrackIndex;
  //   - length
  extern const Var varLongestTrackLength;
  //   - direction
  extern const Var varLongestTrackDirectionX;
  extern const Var varLongestTrackDirectionY;
  extern const Var varLongestTrackDirectionZ;
  extern const Var varLongestTrackDirectionXZ;
  extern const Var varLongestTrackForceDownDirectionX;
  extern const Var varLongestTrackForceDownDirectionY;
  extern const Var varLongestTrackForceDownDirectionZ;
  // - Reco muon


}
