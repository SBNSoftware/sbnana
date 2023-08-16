#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{
  /// \ref SpillCut on valid trigger
  extern const SpillCut kNuMIValidTrigger;

  /// \ref Cut on vertex reconstruced in FV
  extern const Cut kNuMIVertexInFV;

  /// \ref Cut on reco slice not tagged as clearly cosmic
  extern const Cut kNuMINotClearCosmic;

  /// \ref Cut on having muon candidate
  extern const Cut kNuMIHasMuonCandidate;

  /// \ref Cut on having proton candidate
  extern const Cut kNuMIHasProtonCandidate;
  extern const Cut kNuMIProtonCandidateRecoPTreshold;

  /// \ref Cut on having contained primary hadrons
  extern const Cut kNuMIAllPrimaryHadronsContained;

  /// \ref Cut aimed at charged pion rejection
  extern const Cut kNuMINoSecondPrimaryMuonlikeTracks;

  /// \ref Cut aimed at pi0 rejection
  extern const Cut kNuMICutPhotons;

  /// Combined selection \ref Cut for 1muNp0pi with contained+exiting muons
  extern const Cut kNuMISelection_1muNp0pi;

  /// \ref Cut aimed at muon candidate containment, if so desired
  extern const Cut kNuMIMuonCandidateContained;

  /// \ref Cut aimed at reconstruction quality (e.g. split tracks)
  extern const Cut kNuMIRejectSplitMuons;

  /// Signal definitions: Neutrino Neutral Current
  extern const Cut kNuMI_IsSliceNuNC;
  extern const Cut kNuMI_IsSlcNotNu;
  extern const Cut kNuMI_1muNp0piStudy_Signal_NoContainment;
  extern const Cut kNuMI_1muNp0piStudy_OtherNuCC_NoContainment;
  extern const Cut kNuMI_1muNp0piStudy_Signal_NoContainment_ProtonThreshold;
  extern const Cut kNuMI_1muNp0piStudy_OtherNuCC_NoContainment_ProtonThreshold;
}
