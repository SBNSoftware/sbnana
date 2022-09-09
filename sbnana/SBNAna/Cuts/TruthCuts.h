#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{

  extern const Cut kIsAntiNu;
  extern const Cut kIsNu;
  extern const Cut kHasNu;

  extern const Cut kIsNue;
  extern const Cut kIsNumu;
  extern const Cut kIsNutau;

  extern const Cut kIsCC;
  extern const Cut kIsNC;

  extern const Cut kTrueActiveVolumeND;
  extern const Cut kTrueFiducialVolumeND;
  extern const Cut kTrueActiveVolumeFDCryo1;
  extern const Cut kTrueActiveVolumeFDCryo2;

  extern const Cut kTrueFiducialVolumeFDCryo1;
  extern const Cut kTrueFiducialVolumeFDCryo2;

  extern const Cut kVtxDistMagCut;

  extern const Cut kSlcCompletenessCut;

  extern const SpillCut kIsCosmicSpill;
  extern const SpillCut kIsSingleNuSpill;
  extern const SpillCut kIsNueSpill;
  extern const SpillCut kIsNumuSpill;
  extern const SpillCut kIsNutauSpill;

  extern const SpillCut kIsCCSpill;
  extern const SpillCut kIsNCSpill;

  // extern const Cut kIsQE;
  // extern const Cut kIsRes;
  // extern const Cut kIsDIS;
  // extern const Cut kIsCoh;
  // extern const Cut kIsDytmanMEC;


  // //// Select parents
  // extern const Cut kIsPion;
  // extern const Cut kIsKaon;
  // extern const Cut kIsMuon;


  // //// Select neutrinos by their parents type
  // extern const Cut kIsNumuFromKaon;
  // extern const Cut kIsNumuFromPion;
  // extern const Cut kIsNumuFromMuon;

  // extern const Cut kIsNueFromKaon;
  // extern const Cut kIsNueFromPion;
  // extern const Cut kIsNueFromMuon;


  // //// Is the vertex in the detector
  // extern const Cut kIsVtxCont;

}
