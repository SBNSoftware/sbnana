#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Vars/PrimaryUtils.h"

namespace ana
{
  /// bool to determine if object is in fiducial volume
  bool isInFV (double x, double y, double z);

  /// bool to determine if object is in containment volume
  bool isContainedVol (double x, double y, double z);

  /// Utilities for PFParticle loops
  bool IsValidTrkIdx( const caf::SRSliceProxy* slice, const unsigned int idxTrk );
  bool IsTracklikeTrack( const caf::SRSliceProxy* slice, const unsigned int idxTrk );
  bool IsShowerlike( const caf::SRSliceProxy* slice, const unsigned int idxShw );
  bool IsPrimaryPFP( const caf::SRSliceProxy* slice, const unsigned int idxTrk );

  /// \ref Var that is a dummy var that returns 1 for "IsFHC" or 0 for "!IsFHC" for example
  extern const Var kNuMIDummyVar1;
  extern const Var kNuMIDummyVar0;

  /// \ref SpillVar for trigger time (check if the implementation is only comaptible for emulated trigger and fix if so...)
  extern const SpillVar kNuMISpillTriggerTime;

  /// \ref Var for the muon candidate index
  extern const Var kNuMIMuonCandidateIdx;

  /// \ref Var for the proton candidate index
  extern const Var kNuMIProtonCandidateIdx;

  /// \ref MultiVar for the charged pion candidate index
  extern const MultiVar kNuMIChargedPionCandidateIdxs;

  /// \ref MultiVar for the photon candidate indices
  extern const MultiVar kNuMIPhotonCandidateIdxs;

  /// \ref Var for number of showers in event
  extern const Var kNumberRecoShowers;
  
  /// \ref Var for number of tracks in event
  extern const Var kNumberRecoTracks;

  /// kinematic/output variables
  // Muon momentum
  extern const Var kNuMIMuonCandidateRecoP;
  extern const Var kNuMIMuonTrueP;
  // Muon length
  extern const Var kNuMIRecoMuonLength;
  extern const Var kNuMITrueMuonLength;
  // Muon transverse momentum
  extern const Var kNuMIRecoMuonPt;
  extern const Var kNuMITrueMuonPt;
  // Proton momentum
  extern const Var kNuMIProtonCandidateRecoP;
  extern const Var kNuMIProtonTrueP;
  // Proton transverse momentum
  extern const Var kNuMIRecoProtonPt;
  extern const Var kNuMITrueProtonPt;
  // Proton legnth
  extern const Var kNuMIRecoProtonLength;
  extern const Var kNuMITrueProtonLength;
  // Muon angle w.r.t. beam
  extern const Var kNuMIRecoCosThBeam;
  extern const Var kNuMITrueCosThBeam;
  // Muon angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  extern const Var kNuMIRecoCosThVtx;
  extern const Var kNuMITrueCosThVtx;
  // Proton angle w.r.t. beam
  extern const Var kNuMIProtonRecoCosThBeam;
  extern const Var kNuMIProtonTrueCosThBeam;
  // Proton angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  extern const Var kNuMIProtonRecoCosThVtx;
  extern const Var kNuMIProtonTrueCosThVtx;
  // Angle btw muon and proton
  extern const Var kNuMIRecoCosThMuP;
  extern const Var kNuMITrueCosThMuP;

  // TKI variables
  // - delta PT
  extern const Var kNuMIRecodeltaPT;
  extern const Var kNuMITruedeltaPT;
  // - delta PTx
  extern const Var kNuMIRecodeltaPTx;
  extern const Var kNuMITruedeltaPTx;
  // - delta PTy
  extern const Var kNuMIRecodeltaPTy;
  extern const Var kNuMITruedeltaPTy;
  // - delta alphaT
  extern const Var kNuMIRecodeltaalphaT;
  extern const Var kNuMITruedeltaalphaT;
  // - delta phiT
  extern const Var kNuMIRecodeltaphiT;
  extern const Var kNuMITruedeltaphiT;

  // Sideband vars: pi0
  extern const Var kNuMILeadingPhotonCandidateE;
  extern const Var kNuMILeadingPhotonCandidateTrueE;
  extern const Var kNuMISecondaryPhotonCandidateE;
  extern const Var kNuMISecondaryPhotonCandidateTrueE;
  extern const Var kNuMIPhotonCandidatesOpeningAngle;
  extern const Var kNuMILeadingPhotonCandidateLen;
  extern const Var kNuMISecondaryPhotonCandidateLen;
  extern const Var kPi0LeadingPhotonCandidateHitCompletenessBestmatch;
  extern const Var kPi0SecondaryPhotonCandidateHitCompletenessBestmatch;
  extern const Var kPi0LeadingPhotonCandidateEnergyCompletenessBestmatch;
  extern const Var kPi0SecondaryPhotonCandidateEnergyCompletenessBestmatch;
  extern const Var kPi0LeadingPhotonCandidateInFV;
  extern const Var kPi0SecondaryPhotonCandidateInFV;
  extern const Var kPi0LeadingPhotonCandidateBestmatchG4ID;
  extern const Var kPi0SecondaryPhotonCandidateBestmatchG4ID;
  extern const Var kPi0LeadingPhotonCandidateG4ID;
  extern const Var kPi0SecondaryPhotonCandidateG4ID;
  extern const Var kPi0LeadingPhotonCandidateBestplane_Energy;
  extern const Var kPi0SecondaryPhotonCandidateBestplane_Energy;
  extern const Var kPi0LeadingPhotonCandidateBestplane_dEdx;
  extern const Var kPi0SecondaryPhotonCandidateBestplane_dEdx;
  extern const Var kPi0LeadingPhotonCandidateIsContained;
  extern const Var kPi0SecondaryPhotonCandidateIsContained;
  extern const Var kPi0LeadingPhotonCandidateTrackScore;
  extern const Var kPi0SecondaryPhotonCandidateTrackScore;
  extern const Var kPi0LeadingPhotonCandidatePur;
  extern const Var kPi0SecondaryPhotonCandidatePur;
  extern const Var kPi0LeadingPhotonCandidateEff;
  extern const Var kPi0SecondaryPhotonCandidateEff;
  extern const Var kPi0LeadingPhotonCandidateStartX;
  extern const Var kPi0LeadingPhotonCandidateStartY;
  extern const Var kPi0LeadingPhotonCandidateStartZ;
  extern const Var kPi0SecondaryPhotonCandidateStartX;
  extern const Var kPi0SecondaryPhotonCandidateStartY;
  extern const Var kPi0SecondaryPhotonCandidateStartZ;
  extern const Var kPi0LeadingPhotonCandidateTrueStartX;
  extern const Var kPi0LeadingPhotonCandidateTrueStartY;
  extern const Var kPi0LeadingPhotonCandidateTrueStartZ;
  extern const Var kPi0SecondaryPhotonCandidateTrueStartX;
  extern const Var kPi0SecondaryPhotonCandidateTrueStartY;
  extern const Var kPi0SecondaryPhotonCandidateTrueStartZ;
  extern const Var kPi0LeadingPhotonCandidateConversionGap;
  extern const Var kPi0SecondaryPhotonCandidateConversionGap;

  extern const Var kBaryDeltaY;
  extern const Var kBaryDeltaZ;
  extern const Var kBaryRadius;
  extern const Var kBaryFlashFirstHit;

  extern const Var kPi0LeadingPhotonCandidateCosmicDist;
  extern const Var kPi0SecondaryPhotonCandidateCosmicDist;
  extern const Var kPi0LeadingPhotonCandidateShowerDensity;
  extern const Var kPi0SecondaryPhotonCandidateShowerDensity;
}
