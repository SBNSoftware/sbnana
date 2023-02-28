#pragma once

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
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
  extern const SpillVar spillvarNTrack;
  extern const SpillVar spillvarNShower;
  // - CRT Hit
  extern const SpillMultiVar spillvarSideCRTHitPe;
  extern const SpillMultiVar spillvarTopCRTHitPe;
  // - Pos
  extern const SpillMultiVar spillvarEastWestCRTHitPosX;
  extern const SpillMultiVar spillvarFlashPosX;
  // - PMT-CRT matching
  extern const SpillMultiVar spillvarOpFlashTime;
  extern const SpillMultiVar spillvarOpFlashPeakToFirstTime;
  extern const SpillMultiVar spillvarValidOpFlashTime;
  extern const SpillMultiVar spillvarInTimeOpFlashTime; // Valid AND InTime
  extern const SpillMultiVar spillvarInTimeOpFlashPe; // Valid AND InTime, unit : 1e4
  extern const SpillMultiVar spillvarTopCRTHitTime;
  extern const SpillMultiVar spillvarSideCRTHitTime;
  extern const SpillMultiVar spillvarTopCRTPMTTime;
  extern const SpillMultiVar spillvarSideCRTPMTTime;
  extern const SpillMultiVar spillvarCRTPMTTime;
  extern const SpillMultiVar spillvarCRTPMTTimeOfID12;
  extern const SpillMultiVar spillvarCRTPMTMatchingID;
  extern const SpillVar spillvarCRTPMTMatchingEventID;
  // - WW test
  extern const SpillMultiVar spillvarWWTPCTrackEndX;
  extern const SpillMultiVar spillvarWWTPCTrackEndY;
  extern const SpillMultiVar spillvarWWTPCTrackEndZ;

  // Var
  // - GENIE interaction code
  // - https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360
  extern const Var varGENIEIntCode;
  // - Truth interaction
  extern const Var varNeutrinoTruthE;
  extern const Var varNuDirectionX;
  extern const Var varNuDirectionY;
  extern const Var varNuDirectionZ;
  extern const Var varTruthQ2;
  extern const Var varTruthq0_lab;
  extern const Var varTruthmodq_lab;
  extern const Var varTruthW;
  // - Slice var
  extern const Var varCountSlice;
  // - Flash matching
  extern const Var varFMScore;
  extern const Var varFMTime;
  // - NuID
  extern const Var varCRLongestTrackDirY;
  // - Long-enough tracks
  extern const MultiVar varLongTrackDirectionY;
  // - Primary tracks
  extern const MultiVar PrimaryTrackIndices;
  extern const Var NPrimaryTracks;
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
