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
  bool IsProtonLike( const caf::SRSliceProxy* slice, const unsigned int idx );
  bool IsMuonLike( const caf::SRSliceProxy* slice, const unsigned int idx );
  bool IsPionLike( const caf::SRSliceProxy* slice, const unsigned int idx );


  /// \ref Var that is a dummy var that returns 1 for "IsFHC" or 0 for "!IsFHC" for example
  extern const Var kNuMIDummyVar1;
  extern const Var kNuMIDummyVar0;

  /// \ref SpillVar for trigger time (check if the implementation is only comaptible for emulated trigger and fix if so...)
  extern const SpillVar kNuMISpillTriggerTime;
  ///extern const SpillVar kCRTMPMT_FlashMatching;

  /// \ref Var for the muon candidate index
  extern const Var kNuMIMuonCandidateIdx;

  /// \ref Var for the proton candidate index
  extern const Var kNuMIProtonCandidateIdx;

  /// \ref MultiVar for the charged pion candidate index
  extern const MultiVar kNuMIChargedPionCandidateIdxs;

  /// \ref MultiVar for the photon candidate indices
  extern const MultiVar kNuMIPhotonCandidateIdxs;

  /// \ref Var Number of truth matched photons in reconstruction
  extern const Var kNumberTruthMatchRecoPhotons;

  /// \ref Var True opening angle between leading and subleading photon candidates
  extern const Var kNuMITrueCosThPhotonPhoton;
  /// \ref Var True Leading Photon Candidate G4ID
  extern const Var kNuMITrueLeadingPhotonG4ID;
  /// \ref Var True SubLeading Photon Candidate G4ID
  extern const Var kNuMITrueSubLeadingPhotonG4ID;

  /// \ref Var Number of truth matched muons in reconstruction
  extern const Var kNumberTruthMatchRecoMuons;

  /// \ref Var for the leading photon candidate index
  extern const Var kNuMILeadingPhotonCandidateIdx;

  /// \ref Var for the subleading photon candidate index
  extern const Var kNuMISubLeadingPhotonCandidateIdx;

  /// \ref Var for number of showers in the slice
  extern const Var kNumberRecoShowers;
  /// \ref Var for number of pfparticles in the slice
  extern const Var kNumberPFPs;

  /// \ref Var for number of charged pions in slice
  extern const Var kNumberChargedPions;
  
  /// \ref Var for number of tracks in slice
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

  //Reco Nu Vars
  extern const Var kSlcVertexX;
  extern const Var kSlcVertexY;
  extern const Var kSlcVertexZ;

  //Reco Muon Vars
  extern const Var KMuonCandidateRecoStartX;
  extern const Var KMuonCandidateRecoStartY;
  extern const Var KMuonCandidateRecoStartZ;

  //True Muon Vars
  extern const Var kMuonCandidateTrueStartX;
  extern const Var kMuonCandidateTrueStartY;
  extern const Var kMuonCandidateTrueStartZ;
  extern const Var kMuonCandidatePDG;
  extern const Var kProtonCandidatePDG;

  // Sideband vars: pi0
  extern const Var kNuMILeadingPhotonCandidateE;
  extern const Var kNuMILeadingPhotonCandidateTrueE;
  extern const Var kNuMISubLeadingPhotonCandidateE;
  extern const Var kNuMISubLeadingPhotonCandidateTrueE;
  extern const Var kNuMIPhotonCandidatesOpeningAngle;
  extern const Var kNuMILeadingPhotonCandidateLen;
  extern const Var kNuMISubLeadingPhotonCandidateLen;
  extern const Var kPi0LeadingPhotonCandidateHitCompletenessBestmatch;
  extern const Var kPi0SubLeadingPhotonCandidateHitCompletenessBestmatch;
  extern const Var kPi0LeadingPhotonCandidateEnergyCompletenessBestmatch;
  extern const Var kPi0SubLeadingPhotonCandidateEnergyCompletenessBestmatch;
  extern const Var kPi0LeadingPhotonCandidateInFV;
  extern const Var kPi0SubLeadingPhotonCandidateInFV;
  extern const Var kPi0LeadingPhotonCandidateBestmatchG4ID;
  extern const Var kPi0SubLeadingPhotonCandidateBestmatchG4ID;
  extern const Var kPi0LeadingPhotonCandidateG4ID;
  extern const Var kPi0SubLeadingPhotonCandidateG4ID;
  extern const Var kPi0LeadingPhotonCandidateBestplane_Energy;
  extern const Var kPi0SubLeadingPhotonCandidateBestplane_Energy;
  extern const Var kPi0LeadingPhotonCandidateBestplane_dEdx;
  extern const Var kPi0SubLeadingPhotonCandidateBestplane_dEdx;
  extern const Var kPi0LeadingPhotonCandidateIsContained;
  extern const Var kPi0SubLeadingPhotonCandidateIsContained;
  extern const Var kPi0LeadingPhotonCandidateTrackScore;
  extern const Var kPi0SubLeadingPhotonCandidateTrackScore;
  extern const Var kPi0LeadingPhotonCandidatePur;
  extern const Var kPi0SubLeadingPhotonCandidatePur;
  extern const Var kPi0LeadingPhotonCandidateEff;
  extern const Var kPi0SubLeadingPhotonCandidateEff;
  extern const Var kPi0LeadingPhotonCandidateStartX;
  extern const Var kPi0LeadingPhotonCandidateStartY;
  extern const Var kPi0LeadingPhotonCandidateStartZ;
  extern const Var kPi0SubLeadingPhotonCandidateStartX;
  extern const Var kPi0SubLeadingPhotonCandidateStartY;
  extern const Var kPi0SubLeadingPhotonCandidateStartZ;
  extern const Var kPi0LeadingPhotonCandidateTrueStartX;
  extern const Var kPi0LeadingPhotonCandidateTrueStartY;
  extern const Var kPi0LeadingPhotonCandidateTrueStartZ;
  extern const Var kPi0SubLeadingPhotonCandidateTrueStartX;
  extern const Var kPi0SubLeadingPhotonCandidateTrueStartY;
  extern const Var kPi0SubLeadingPhotonCandidateTrueStartZ;
  extern const Var kPi0LeadingPhotonCandidateConversionGap;
  extern const Var kPi0SubLeadingPhotonCandidateConversionGap;

  extern const Var kBaryDeltaY;
  extern const Var kBaryDeltaZ;
  extern const Var kBaryRadius;
  extern const Var kBaryFlashFirstHit;

  extern const Var kPi0LeadingPhotonCandidateCosmicDist;
  extern const Var kPi0SubLeadingPhotonCandidateCosmicDist;
  extern const Var kPi0LeadingPhotonCandidateShowerDensity;
  extern const Var kPi0SubLeadingPhotonCandidateShowerDensity;

  extern const Var kPi0LeadingPhotonCandidateShowerGenPX;
  extern const Var kPi0LeadingPhotonCandidateShowerGenPY;
  extern const Var kPi0LeadingPhotonCandidateShowerGenPZ;
  extern const Var kPi0SubLeadingPhotonCandidateShowerGenPX;
  extern const Var kPi0SubLeadingPhotonCandidateShowerGenPY;
  extern const Var kPi0SubLeadingPhotonCandidateShowerGenPZ;
  extern const Var kPi0LeadingPhotonCandidateTrueCryostat;
  extern const Var kPi0SubLeadingPhotonCandidateTrueCryostat;
  extern const Var kPi0LeadingPhotonCandidateTrueLength;
  extern const Var kPi0SubLeadingPhotonCandidateTrueLength;
  extern const Var kPi0LeadingPhotonCandidatePDG;
  extern const Var kPi0SubLeadingPhotonCandidatePDG;
  extern const Var kPi0LeadingPhotonCandidateProtonChi2;
  extern const Var kPi0SubLeadingPhotonCandidateProtonChi2;
  extern const Var kPi0LeadingPhotonCandidateMuonChi2;
  extern const Var kPi0SubLeadingPhotonCandidateMuonChi2;
  extern const Var kPi0LeadingPhotonCandidatePionChi2;
  extern const Var kPi0SubLeadingPhotonCandidatePionChi2;
  extern const Var kPi0LeadingPhotonCandidateMuonPionChi2Diff;
  extern const Var kPi0SubLeadingPhotonCandidateMuonPionChi2Diff;
  extern const Var kPi0LeadingPhotonCandidateNHits;
  extern const Var kPi0SubLeadingPhotonCandidateNHits;
  extern const Var kPi0LeadingPhotonCandidateSqrtEDen;
  extern const Var kPi0SubLeadingPhotonCandidateSqrtEDen;
  extern const Var kPi0LeadingPhotonCandidateHitDen;
  extern const Var kPi0SubLeadingPhotonCandidateHitDen;

  extern const Var kIsClearCosmic;

}
