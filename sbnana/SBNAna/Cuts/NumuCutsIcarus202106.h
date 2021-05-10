#pragma once

// Definition of the generic Cut object.
#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{
  // Cuts for identifying slice category.
  extern const Cut kIsNuSlice;

  extern const Cut kIsCosmic;

  extern const Cut kIsNuMuCC;

  extern const Cut kIsNuOther;

  // "Raw" cuts for the selection.
  extern const Cut kCryo0;

  extern const Cut kTFiducial;

  extern const Cut kRFiducial;

  extern const Cut kNotClearCosmic;

  extern const Cut kNuScore;

  extern const Cut kFMScore;

  extern const Cut kPTrack;

  extern const Cut kPTrackContained;
  
  extern const Cut kPTrackExiting;

  // "Successive" cuts for the selection.
  extern const Cut kNuMuCC_Cryo0;

  extern const Cut kCosmic_Cryo0;

  extern const Cut kNuOther_Cryo0;

  extern const Cut kNuMuCC_TFiducial;

  extern const Cut kNuMuCC_RFiducial;

  extern const Cut kCosmic_TFiducial;

  extern const Cut kCosmic_RFiducial;

  extern const Cut kNuOther_TFiducial;

  extern const Cut kNuOther_RFiducial;

  extern const Cut kNuMuCC_ClearCos;

  extern const Cut kCosmic_ClearCos;

  extern const Cut kNuOther_ClearCos;

  extern const Cut kNuMuCC_NuScore;

  extern const Cut kCosmic_NuScore;

  extern const Cut kNuOther_NuScore;

  extern const Cut kNuMuCC_FMScore;

  extern const Cut kCosmic_FMScore;

  extern const Cut kNuOther_FMScore;

  extern const Cut kNuMuCC_PTrack;

  extern const Cut kCosmic_PTrack;

  extern const Cut kNuOther_PTrack;

  // The full selection cut using only reco information.
  extern const Cut kNuMuCC_FullSelection;

  // Spill level cuts.
  extern const SpillCut kLongTrack;
}
