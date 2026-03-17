//! ////////////////////////////////////////////////////////////////////////////
//! @file: nuESelection_Cuts.h                                                //
//! @authors: Jacob Smith (smithja)                                           //
//! @email: jacob.a.smith@stonybrook.edu                                      //
//! Last edited: November 15th, 2025                                          //
//!                                                                           //
//! @brief Header file to define cuts for the ICARUS electron neutrino        //
//! selection. Accompanying information about plot formatting (e.g. plot      //
//! colors) can be found in nuESelectionHelper_icarus.h.                      //
//! ////////////////////////////////////////////////////////////////////////////
#pragma once

#include <cmath>  //! to use std::isnan
#include <iostream>    //! to use std::cout and std::cin, etc.

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h" //! to work with CAF data

//! CAFAna objects used in this file:
#include "sbnana/CAFAna/Core/Cut.h"

//! SBNAna objects used in this file:
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NueCuts.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/NueVars.h"

//! Imports user-defined general tools developed over a research career:
#include "auxTools.h"

namespace ana {

//! ////////////////////////////////////////////////////////////////////////////
//! Nominal selection cut values found by manual tuning:
struct NominalCuts {
  static constexpr double showerEnergy_lb = 0.3;   //! [GeV]
  static constexpr double showerdEdx_ub = 3.0;     //! [MeV/cm]
  static constexpr double showerConvGap_ub = 3.25; //! [cm]
  static constexpr double showerDensity_lb = 4.5;  //! [MeV/cm]
  static constexpr double nuETrackLen_ub = 110.0;  //! [cm]
  static constexpr double flashMatch_lb = 0.0;     //! [dimensionless]
  static constexpr double flashMatch_ub = 6.0;     //! [dimensionless]
  static constexpr double barycenter_lb = 0.0;     //! [cm]
  static constexpr double barycenter_ub = 30.0;    //! [cm]
  static constexpr double trackScore_ub = 0.5;     //! [dimensionless]
};
//! ////////////////////////////////////////////////////////////////////////////

//! ////////////////////////////////////////////////////////////////////////////
//! Cut and SpillCut Definitions (inline for header-only usage):
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
//! Aliased Cuts from SBNAna (for naming conventions):
//! Further details on how/where these Cuts are created can be found in
//! NueCuts.cxx and NueCuts.h.
//!
//! 1. kRecoShower is a rough quality cut ensuring non-zero best-plane energy,
//!    conversion gap, and dE/dx.
//! 2. kRecoShowerFD comes from NueCuts.h. @note this doesn't include any
//!    cut on a raw track score variable:
//!        const Cut kRecoShowerFD = kRecoShower && kNueNumShowersCut \
//!                                 && kShowerdEdxCut && kShowerConvGapCut \
//!                                 && kNueTrackLenCut && kShowerDensityCut \
//!                                 && kShowerEnergyCut;
//! 3. kNueContainedFD ensures the shower is contained in fiducial volume.
//! 4. kNueNumShowersCut ensures there's exactly one shower in the slice.
//! 5. kNueTrackLenCut ensures longest track is less than 110 cm.

const Cut kRecoShowerCut = kRecoShower; //! add "Cut" suffix for clarity
const Cut kRecoShowerFDCut = kRecoShowerFD; //! add "Cut" suffix for clarity
const Cut kNuEContainedFDCut = kNueContainedFD; //! add "Cut" suffix for clarity
                                                //! and change "Nue" to "NuE"
const Cut kNuENumShowersCut = kNueNumShowersCut; //! change "Nue" to "NuE"
const Cut kNuETrackLenCut = kNueTrackLenCut; //! change "Nue" to "NuE"

//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
//! Individual Cuts and SpillCuts:
const Cut kBarycenterFMDeltaZCut([](const caf::SRSliceProxy *slc) {
  return !std::isnan(slc->barycenterFM.deltaZ_Trigger) \
          && slc->barycenterFM.deltaZ_Trigger >= 0 \
          && slc->barycenterFM.deltaZ_Trigger < 100;
});

const SpillCut kContainedSpillCut([](const caf::SRSpillProxy* sr) {
  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kNuEContainedFDCut(&slc))
      ++counter;
  }

  return counter > 0;
});

const SpillCut kFlashMatchSpillCut([](const caf::SRSpillProxy* sr) {
  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kSlcFlashMatchCut(&slc))
      ++counter;
  }

  return counter > 0;
});


const SpillCut kRecoShowerFDSpillCut([](const caf::SRSpillProxy* sr) {
  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kRecoShowerFDCut(&slc))
      ++counter;
  }

  return counter > 0;
});

//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
//! CUSTOM Cut factory functions (return Cuts with user-specified thresholds):
inline Cut makeCUSTOMShowerEnergyCut(double lb) {
  return Cut([lb](const caf::SRSliceProxy *slc) {
    const int largestShwIdx(kLargestRecoShowerIdx(slc));
    if (largestShwIdx == -1) return false;
    return (slc->reco.pfp[largestShwIdx].shw.bestplane_energy > lb);
  });
}

inline Cut makeCUSTOMShowerdEdxCut(double ub) {
  return Cut([ub](const caf::SRSliceProxy *slc) {
    const int largestShwIdx(kLargestRecoShowerIdx(slc));
    if (largestShwIdx == -1) return false;
    return (slc->reco.pfp[largestShwIdx].shw.bestplane_dEdx < ub);
  });
}

inline Cut makeCUSTOMShowerConvGapCut(double ub) {
  return Cut([ub](const caf::SRSliceProxy *slc) {
    const int largestShwIdx(kLargestRecoShowerIdx(slc));
    if (largestShwIdx == -1) return false;
    return (slc->reco.pfp[largestShwIdx].shw.conversion_gap < ub);
  });
}

inline Cut makeCUSTOMShowerDensityCut( double lb) {
  return Cut([lb](const caf::SRSliceProxy* slc) {
    const int largestShwIdx(kLargestRecoShowerIdx(slc));
    if ( largestShwIdx==-1 ) return false;
    return (slc->reco.pfp[largestShwIdx].shw.density > lb);
  });
}

inline Cut makeCUSTOMNuETrackLenCut( double ub) {
  return Cut([ub](const caf::SRSliceProxy *slc) {
    return kLongestTrackLength(slc) < ub;
  });
}

inline Cut makeCUSTOMSlcFlashMatchCut(double lb, double ub) {
  return Cut([lb, ub](const caf::SRSliceProxy *slc) {
    return (kSlcHasFlashMatch(slc) && slc->fmatch.score > lb && slc->fmatch.score < ub);
  });
}

inline Cut makeCUSTOMBarycenterFMDeltaZCut(double lb_eq, double ub) {
  return Cut([lb_eq, ub](const caf::SRSliceProxy *slc) {
    return !std::isnan(slc->barycenterFM.deltaZ_Trigger) &&
            slc->barycenterFM.deltaZ_Trigger >= lb_eq &&
            slc->barycenterFM.deltaZ_Trigger < ub;
  });
}

inline Cut makeCUSTOMTrackScoreCut(double ub) {
  return Cut([ub](const caf::SRSliceProxy *slc) {
    for( unsigned int i=0; i < slc->reco.npfp; i++) {
      if( slc->reco.pfp[i].trackScore > ub) return false;
    }
    return true;
  });
}

//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
//! Truth-level Cuts:
const Cut kTrueNuECC      = kIsNue  && !kIsNC; //! charged-current electron neutrino
const Cut kTrueNuMuCC     = kIsNumu && !kIsNC; //! charged-current muon neutrino
const Cut kTrueNC         = kIsNC;             //! neutral current interactions
const Cut kTrueCosmic     = !kHasNu;           //! no neutrinos? then it's a cosmic!

//! kTrueFiducialVolumeFD cuts: kHasNu && vertex in specified volume
const Cut kTrueFVFDCut = kTrueFiducialVolumeFDCryo1 || kTrueFiducialVolumeFDCryo2;

//! @note: The following SpillCuts all return false if the sr is a cosmic or
//!        single-neutrino spill:
const SpillCut kTrueNuECCSpillCut    = kIsNueSpill && !kIsNCSpill;  //! std::abs( sr->mc.nu[0].pdg) == 12 and !sr->mc.nu[0].isnc
const SpillCut kTrueNuMuCCSpillCut   = kIsNumuSpill && !kIsNCSpill; //! std::abs( sr->mc.nu[0].pdg) == 14 and !sr->mc.nu[0].isnc
const SpillCut kTrueNCSpillCut       = kIsNCSpill;                  //! sr->mc.nu[0].isnc
const SpillCut kTrueCosmicSpillCut   = kIsCosmicSpill;              //! sr->mc.nnu == 0 ub) {

//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
//! Reconstruction-level Cuts:

//! @note: kFiducialVolumeFD cuts require vertex in specified volume
const Cut kRecoFVFDCut = kFiducialVolumeFDCryo1 || kFiducialVolumeFDCryo2;

//! @note: We define a "true signal" for electron neutrino interactions.
const Cut kSignal = kRecoFVFDCut && kNuEContainedFDCut && kShowerEnergyCut;

//! @note: kSlcIsRecoNu = !slc->is_clear_cosmic is the reconstruction-level
//! quantitiy, whereas kTrueNuECC, kTrueNuMuCC, kTrueNC, and kTrueCosmic are
//! truth-level quantities.
const Cut kRecoNuECC = kSlcIsRecoNu && kTrueNuECC;
const Cut kRecoNuECCSignal = kSlcIsRecoNu && kTrueNuECC && kSignal;
const Cut kRecoNuMuCC = kSlcIsRecoNu && kTrueNuMuCC;
const Cut kRecoNC = kSlcIsRecoNu && kTrueNC;
const Cut kRecoCosmic = kSlcIsRecoNu && kTrueCosmic;

//! Shorthand for non-reco shower cut used in some of N-1 definitions
const Cut kNonRecoShowerFD = kNuEContainedFDCut && kSlcFlashMatchCut && kBarycenterFMDeltaZCut;

//! @note: N-1 cuts are manually defined since we have an overall reconstruction
//! cut that is comprised of other cuts listed in individual_cuts_slice.
const Cut kN1ContainedCut          =                       kSlcFlashMatchCut && kBarycenterFMDeltaZCut && kRecoShowerFDCut;
const Cut kN1FlashMatchCut         = kNuEContainedFDCut                      && kBarycenterFMDeltaZCut && kRecoShowerFDCut;
const Cut kN1BarycenterFMDeltaZCut = kNuEContainedFDCut && kSlcFlashMatchCut                           && kRecoShowerFDCut;
const Cut kN1RecoShowerFDCut       = kNuEContainedFDCut && kSlcFlashMatchCut && kBarycenterFMDeltaZCut;

const Cut kN1RecoShowerCut    = kNonRecoShowerFD                   && kShowerEnergyCut && kShowerdEdxCut && kShowerConvGapCut && kShowerDensityCut && kNuETrackLenCut && kNuENumShowersCut;
const Cut kN1ShowerEnergyCut  = kNonRecoShowerFD && kRecoShowerCut                     && kShowerdEdxCut && kShowerConvGapCut && kShowerDensityCut && kNuETrackLenCut && kNuENumShowersCut;
const Cut kN1ShowerdEdxCut    = kNonRecoShowerFD && kRecoShowerCut && kShowerEnergyCut                   && kShowerConvGapCut && kShowerDensityCut && kNuETrackLenCut && kNuENumShowersCut;
const Cut kN1ShowerConvGapCut = kNonRecoShowerFD && kRecoShowerCut && kShowerEnergyCut && kShowerdEdxCut                      && kShowerDensityCut && kNuETrackLenCut && kNuENumShowersCut;
const Cut kN1ShowerDensityCut = kNonRecoShowerFD && kRecoShowerCut && kShowerEnergyCut && kShowerdEdxCut && kShowerConvGapCut                      && kNuETrackLenCut && kNuENumShowersCut;
const Cut kN1TrackLenCut      = kNonRecoShowerFD && kRecoShowerCut && kShowerEnergyCut && kShowerdEdxCut && kShowerConvGapCut && kShowerDensityCut                    && kNuENumShowersCut;
const Cut kN1NumShowersCut    = kNonRecoShowerFD && kRecoShowerCut && kShowerEnergyCut && kShowerdEdxCut && kShowerConvGapCut && kShowerDensityCut && kNuETrackLenCut;

const Cut kFullCut = kNuEContainedFDCut \
                      && kSlcFlashMatchCut \
                      && kBarycenterFMDeltaZCut \
                      && kRecoShowerFDCut;
//! ////////////////////////////////////////////////////////////////////////////
} //! end namespace ana