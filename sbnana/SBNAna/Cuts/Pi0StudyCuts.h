#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{
  /// \ref SpillCut on valid trigger
  extern const SpillCut kNuMIValidTrigger;
  extern const SpillCut kBNBValidTrigger;

  /// \ref Cut on vertex reconstruced in FV
  extern const Cut kNuMIVertexInFV;

  /// \ref Cut on vertex reconstruced in contained Volume
  extern const Cut kNuMIVertexIsContained;

  /// \ref Cut on reco slice not tagged as clearly cosmic
  extern const Cut kNuMINotClearCosmic;

  /// \ref Cut on having muon candidate
  extern const Cut kNuMIHasMuonCandidate;

  /// \ref Cut on having proton candidate
  extern const Cut kNuMIHasProtonCandidate;
  //extern const Cut kNuMIProtonCandidateRecoPTreshold;

  /// \ref Cut on the barycenter flash
  extern const Cut kBaryDeltaZTCut;
  extern const Cut kBaryRadiusTCut;
  extern const Cut kBaryFlashFHitCut;

  /// \ref Cut on having contained primary hadrons
  extern const Cut kNuMIAllPrimaryHadronsContained;

  /// \ref Cut aimed at charged pion rejection
  extern const Cut kNuMINoSecondPrimaryMuonlikeTracks;

  /// \ref Cut aimed at pi0 rejection
  extern const Cut kNuMICutPhotons;

  ///Base selection of 1muXPi0 for the main 1mu1Pi0 selection
  extern const Cut kNuMISelection_1muXpi0_Base;

  /// added cosmics and charged pion filters
  extern const Cut kNuMISelection_1muXpi0;
  
  extern const Cut kNuMISelection_1muXpi0_CosmicFilter;
  extern const Cut kNuMISelection_1muXpi0_ChargedPionFilter;


  /// Base selection of 1muNp from which the main selection and current sidebands branch
  extern const Cut kNuMISelection_1muNp_Base;

  /// Combined selection \ref Cut for 1muNp0pi with contained+exiting muons, without additional shower cut being used
  extern const Cut kNuMISelection_1muNp0pi_WithoutShowerCut;

  /// Combined selection \ref Cut for 1muNp0pi with contained+exiting muons
  extern const Cut kNuMISelection_1muNp0pi;

  /// \ref Cut aimed at muon candidate containment, if so desired
  extern const Cut kNuMIMuonCandidateContained;

  /// \ref Cut aimed at reconstruction quality (e.g. split tracks)
  extern const Cut kNuMIRejectSplitMuons;

  /// \ref Cut aimed at having TWO photons for better pi0 selection
  extern const Cut kNuMIHasTwoPhotons;

  /// \ref Cut pion sideband
  extern const Cut kNuMIChargedPionSideBand;
  extern const Cut kNuMINeutralPionSideBand;
  extern const Cut kNuMINeutralPion2phSideBand;

  /// \ref CutType; 1=Signal, 2=pi+- sideband, 3=pi0 sideband (0=other)
  extern const Var kNuMICutType;
  extern const Var kNuMICutTypeWithoutShowerCut;

  /// \ref Signal definitions: Neutrino Neutral Current
  extern const Cut kNuMI_IsSliceNuNC;
  /// \ref Not nu matched: i.e. cosmic, or noise, or not well-matched to an interaction
  extern const Cut kNuMI_IsSlcNotNu;
  /// \ref Check 1muNp0pi using vector of primaries
  bool Is1muNp0pi(const caf::Proxy<caf::SRTrueInteraction>& true_int, bool ApplyProtonPCut);
  bool Is1mu1Pi0X(const caf::Proxy<caf::SRTrueInteraction>& true_int); //Pi0Signal
  inline bool Is1muNp0piWithProtonPcut(const caf::Proxy<caf::SRTrueInteraction>& true_int){ return Is1muNp0pi(true_int, true); }
  extern const Cut kNuMI_1muNp0piStudy_Signal_NoContainment;
  extern const Cut kNuMI_1mu1Pi0XStudy_Signal_NoContainment;
  extern const Cut kNuMI_1muNp0piStudy_OtherNuCC_NoContainment;
  extern const Cut kNuMI_1muNp0piStudy_Signal_NoContainment_ProtonThreshold;
  extern const Cut kNuMI_1muNp0piStudy_OtherNuCC_NoContainment_ProtonThreshold;

  /// \ref Var for slice type (signal, other NuCC, NuNC, NotNu)
  extern const Var kNuMISliceSignalType;

}
