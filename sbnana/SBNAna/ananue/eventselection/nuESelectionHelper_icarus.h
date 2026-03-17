//! ////////////////////////////////////////////////////////////////////////////
//! @file: nuESelectionHelper_icarus.h                                        //
//! @authors: Jacob Smith (smithja)                                           //
//! @email: jacob.a.smith@stonybrook.edu                                      //
//! Last edited: November 9th, 2025                                           //
//!                                                                           //
//! @details: Information about electron neutrino selection                   //
//! and associated plot formatting (e.g. plot colors).                        //
//! ////////////////////////////////////////////////////////////////////////////
#pragma once

#include <vector>      //! to use std::vector objects
#include <string>      //! to use std::string objects
#include <iostream>    //! to use std::cout and std::cin, etc.
#include <limits>      //! to use std::numeric_limits

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h" //! to work with CAF data

//! CAFAna objects used in this file:
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Var.h"

//! SBNAna objects used in this file:
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Vars/Binnings.h"

#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"

#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Cuts/NueCuts.h"

//! Additional Cuts specific to ICARUS electron neutrino selection:
#include "nuESelection_Cuts.h"

//! Imports user-defined general tools developed over a research career:
#include "auxTools.h"

using namespace ana;


//! Declare canvas resolutions:
Int_t cXRes = 700; //! resolution of width of output canvas
Int_t cYRes = 500; //! resolution of height of output canvas




//! Alias interaction type colors:
Color_t colorEff     = kGreen + 2;   //! efficiency
Color_t colorPur     = kOrange + 2;  //! purity
Color_t colorTotal   = kBlack;       //! anything with a reconstructed neutrino
Color_t colorNuE     = kBlue - 7;    //! electron neutrino
Color_t colorNuMu    = kGreen + 1;   //! muon neutrino
Color_t colorNC      = kMagenta;     //! neutral current interactions
Color_t colorCosmics = kGray + 2;    //! cosmics
Color_t colorOther   = kOrange + 8;  //! anything that's not NuECC, NuMuCC, NC,
                                     //! or cosmic

//! Alias selection cut colors:
Color_t colorNoCut                 = kBlack;      //! kNoCut
Color_t colorContainedCut          = kBlue;       //! kNuEContainedFDCut
Color_t colorFlashMatchCut         = kOrange + 8; //! kSlcFlashMatchCut and kFlashMatchSpillCut
Color_t colorRecoShowerFDCut       = kViolet + 7; //! kRecoShowerFDCut and kRecoShowerFDSpillCut
Color_t colorBarycenterFMDeltaZCut = kGreen + 3;  //! kBarycenterFMDeltaZCut

//! @note: These Cuts are constituent Cuts of RecoShowerFDCut.
//!        So we have them be similar in color to colorRecoShowerFDCut.
Color_t colorRecoShowerCut     = kViolet - 1; //! kRecoShowerCut
Color_t colorNuENumShowersCut  = kViolet - 2; //! kNuENumShowersCut
Color_t colorShowerdEdxCut     = kViolet - 3; //! kShowerdEdxCut
Color_t colorShowerConvGapCut  = kViolet - 4; //! kShowerConvGapCut
Color_t colorNuETrackLenCut    = kViolet - 5; //! kNueTrackLenCut
Color_t colorShowerDensityCut  = kViolet - 6; //! kShowerDensityCut
Color_t colorShowerEnergyCut   = kViolet - 7; //! kShowerEnergyCut

Color_t colorCRTHitVetoFDSpillCut = kViolet - 8; //! CRTHitVeto (SpillCut only!)
Color_t colorFullCut              = kRed;        //! all selection cuts


//! Alias line styles:
Style_t lineNuE        = kSolid;
Style_t lineTrueSigNuE = kDashed;
Style_t lineNuMu       = kSolid;
Style_t lineNC         = kDashed;
Style_t lineCosmics    = kSolid;
Style_t lineOther      = kDotted;


//! ////////////////////////////////////////////////////////////////////////////
//! Binning Definitions: units for corresponding variables are given in
//! Variables section
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
const Binning kEnergyBinning        = Binning::Simple( 120,   0.0,    3.0);
const Binning kDedxBinning          = Binning::Simple( 40,    0.0,    10);
const Binning kGapBinning           = Binning::Simple( 40,    0.0,    10);
const Binning kDensityBinning       = Binning::Simple( 40,    0.0,    10);
const Binning kOpenAngleBinning     = Binning::Simple( 60,    0.0,    180);
const Binning kLengthBinning        = Binning::Simple( 40,    0.0,    200);
const Binning kShowerEnergyBinning  = Binning::Simple( 120,   0.0,    3.0);
const Binning kBarycenterFMBinning  = Binning::Simple( 41,   -5.0,    200.0);
const Binning kFMScoreBinning       = Binning::Simple( 30,    0.0,    15.0);

//! CRT Variables:
const Binning kPEBinning            = Binning::Simple( 60,   0.0,    600);
const Binning kTimeBinning          = Binning::Simple( 155, -1550.0, 1550.0);
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Here we define the vectors of variables, selection cuts, and interaction
//! types that we will investigate as part of this electron neutrino selection
//! analysis:

/** @struct PlotDefSlice
 *  @brief Information used in slice-level variables.
 */
struct PlotDefSlice
{
  std::string label;
  std::string unit;
  Binning bins;
  Var var;
};

/** @struct PlotDefSpill
 *  @brief Information used in spill-level variables.
 */
struct PlotDefSpill
{
  std::string label;
  std::string unit;
  Binning bins;
  SpillVar var;
};

/** @struct PlotDefMultiVar
 *  @brief Information used in MultiVar-type variables.
 */
struct PlotDefMultiVar
{
  std::string label;
  std::string unit;
  Binning bins;
  SpillMultiVar var;
};

/** @struct SelDefSlice
 *  @brief Information used in slice-level cuts.
 */
struct SelDefSlice
{
  std::string label;
  Cut cut;
  int color;
  std::string desc; //! short description
};

/** @struct SelDefSpill
 *  @brief Information used in spill-level cuts.
 */
struct SelDefSpill
{
  std::string label;
  SpillCut cut;
  int color;
  std::string desc; //! short dscription
};



//! ////////////////////////////////////////////////////////////////////////////
//! Variables:
//! @note: empty second parameter implies no units/dimensions
std::vector<PlotDefSlice> vars_slice = {
  {"NumberOfSlices",            "",           Binning::Simple(3,0,3), kCounting},
  {"RecoShowerOpeningAngle",    "deg",        kOpenAngleBinning,      kRecoShower_OpenAngle},
  {"RecoShowerStartX",          "cm",         kPositionXFDBinning,    kRecoShower_StartX},
  {"RecoShowerStartY",          "cm",         kPositionYFDBinning,    kRecoShower_StartY},
  {"RecoShowerStartZ",          "cm",         kPositionZFDBinning,    kRecoShower_StartZ},
  {"RecoShowerEndX",            "cm",         kPositionXFDBinning,    kRecoShower_EndX},
  {"RecoShowerEndY",            "cm",         kPositionYFDBinning,    kRecoShower_EndY},
  {"RecoShowerEndZ",            "cm",         kPositionZFDBinning,    kRecoShower_EndZ},
  {"SliceVertexX",              "cm",         kPositionXFDBinning,    kSlcVtxX},
  {"SliceVertexY",              "cm",         kPositionYFDBinning,    kSlcVtxY},
  {"SliceVertexZ",              "cm",         kPositionZFDBinning,    kSlcVtxZ},
  {"RecoShowerConversionGap",   "cm",         kGapBinning,            kRecoShower_ConversionGap},
  {"RecoShowerBestPlane_dEdx",  "MeV_per_cm", kDedxBinning,           kRecoShower_BestdEdx},
  {"RecoShowerBestPlaneEnergy", "GeV",        kEnergyBinning,         kRecoShower_BestEnergy},
  {"RecoShowerDensity",         "MeV_per_cm", kDensityBinning,        kRecoShower_Density},
  {"RecoShowerEnergy",          "GeV",        kEnergyBinning,         kRecoShower_Energy},
  {"RecoShowerLength",          "cm",         kLengthBinning,         kRecoShower_Length},
  {"LongestTrackLength",        "cm",         kLengthBinning,         kLongestTrackLength},
  {"BarycenterDeltaZTrigger",   "cm",         kBarycenterFMBinning,   kBarycenterFMDeltaZ},
  {"FlashScore_SLICE",          "",           kFMScoreBinning,        kSlcFlashScore},
  {"TrueEnergy",                "GeV",        kLowEnergyGeVBinning,   kTruthEnergy}
};

std::vector<PlotDefSpill> vars_spill = {
  {"NumberOfSpills",            "",    Binning::Simple(3,0,3), kSpillCounting},
  {"TrueNeutrinoEnergyInSpill", "GeV", kLowEnergyGeVBinning,   kTruthNuEnergy},
  {"TrueLeptonEnergyInSpill",   "GeV", kLowEnergyGeVBinning,   kTruthLeptonEnergy}
};

//! @note: There are no crtvars_slice entries since CRT variables are all at
//!        the spill level.
std::vector<PlotDefMultiVar> crtvars_spill = {
  {"CRTHitX",    "cm", kCRTXFDBinning, kCRTHitX},
  {"CRTHitY",    "cm", kCRTYFDBinning, kCRTHitY},
  {"CRTHitZ",    "cm", kCRTZFDBinning, kCRTHitZ},
  {"CRTHitTime", "us", kTimeBinning,   kCRTHitTimeFD},
  {"CRTHitPE",   "PE", kPEBinning,     kCRTHitPE}
};
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Selection Cuts:
std::vector<SelDefSlice> individual_cuts_slice = {
  {"kNoCut", kNoCut, colorNoCut, "NO SELECTION CUTS APPLIED."},
  //! Low-level reconstruction Cuts:
  {"kRecoShower",    kRecoShowerCut,    colorRecoShowerCut,    "Largest shower in this slice has positive energy, dEdx, and conversion gap."},
  {"kShowerEnergy",  kShowerEnergyCut,  colorShowerEnergyCut,  "Largest shower in this slice has an energy greater than 200 MeV."},
  {"kShower_dEdx",   kShowerdEdxCut,    colorShowerdEdxCut,    "Largest shower in this slice has dEdx less than 3.625 MeV per cm."},
  {"kShowerConvGap", kShowerConvGapCut, colorShowerConvGapCut, "Largest shower in this slice has a conversion gap less than 3.25 cm."},
  {"kShowerDen",     kShowerDensityCut, colorShowerDensityCut, "Largest shower in this slice has a density greater than 4.5 MeV per cm."},
  {"kNuETrackLen",   kNueTrackLenCut,   colorNuETrackLenCut,   "Longest track in this slice is less than 110 cm."},
  //! Containment and Fiducial Volume Cuts:
  {"kNuEContainedFD", kNuEContainedFDCut, colorContainedCut, "Largest shower in this slice is fully contained."},
  //! Particle ID and High-level reconstruction Cuts:
  {"kSlcFlashMatch",      kSlcFlashMatchCut,      colorFlashMatchCut,         "Slice FlashMatch score within (0, 6)."},
  {"kBarycenterFMDeltaZ", kBarycenterFMDeltaZCut, colorBarycenterFMDeltaZCut, "Slice barycenterFM.deltaZ_Trigger within [0, 100)."},
  {"kNuENumShowers",      kNuENumShowersCut,      colorNuENumShowersCut,      "Reconstructed PFP in this slice has exactly one shower with energy greater than 200 MeV."},
  //! Conglomerate Cuts:
  {"kRecoShowerFDCut", kRecoShowerFDCut, colorRecoShowerFDCut, "The following Cuts: kRecoShower, kShowerEnergy, kShowerdEdx, kShowerConvGap, " \
                                                                       " kShowerDensity, kNueTrackLen, kNueNumShowers."},
  {"FullSelection", kFullCut, colorFullCut, "ALL SELECTION CUTS APPLIED."}
};
const size_t nCuts = individual_cuts_slice.size();

std::vector<SelDefSlice> nMinus1_cuts_slice = {
  {"kN1RecoShowerCut",       kN1RecoShowerCut,         colorRecoShowerCut,         "Full Selection minus minimum thresholds on largest shower energy, dEdx, and conversion gap."},
  {"kN1ShowerEnergy",        kN1ShowerEnergyCut,       colorShowerEnergyCut,       "Full Selection minus maximum threshold on largest shower energy."},
  {"kN1Shower_dEdx",         kN1ShowerdEdxCut,         colorShowerdEdxCut,         "Full Selection minus maximum threshold on largest shower dEdx."},
  {"kN1ShowerConvGap",       kN1ShowerConvGapCut,      colorShowerConvGapCut,      "Full Selection minus maximum threshold on largest shower conversion gap."},
  {"kN1ShowerDen",           kN1ShowerDensityCut,      colorShowerDensityCut,      "Full Selection minus maximum threshold on largest shower density."},
  {"kN1NuETrackLen",         kN1TrackLenCut,           colorNuETrackLenCut,        "Full Selection minus maximum threshold on longest track in slice."},
  {"kN1Contained",           kN1ContainedCut,          colorContainedCut,          "Full Selection minus shower is contained in FV."},
  {"kN1SlcFlashMatch",       kN1FlashMatchCut,         colorFlashMatchCut,         "Full Selection minus slice FlashMatch score within (0, 6)."},
  {"kN1BarycenterFMDeltaZ",  kN1BarycenterFMDeltaZCut, colorBarycenterFMDeltaZCut, "Full Selection minus barycenterFM.deltaZ_Trigger within [0, 100)."},
  {"kN1NuENumShowers",       kN1NumShowersCut,         colorNuENumShowersCut,      "Full Selection minus slice has exactly one shower with energy greater than 0.2 MeV."},
  {"kN1RecoShowerFDCut",     kN1RecoShowerFDCut,       colorRecoShowerFDCut,       "Functionally, the containment, FlashMatch score, and barycenterFM.deltaZ_Trigger Cuts."}
};

//! Add regular and N-1 cuts into a conglomerate selections vector:
std::vector<SelDefSlice> sels_slice = [] {
  std::vector<SelDefSlice> sels_slice;
  // Add all individual cuts
  sels_slice.insert(sels_slice.end(),
                    individual_cuts_slice.begin(),
                    individual_cuts_slice.end());

  // Add all N–1 cuts
  sels_slice.insert(sels_slice.end(),
                    nMinus1_cuts_slice.begin(),
                    nMinus1_cuts_slice.end());

  return sels_slice;
}();

std::vector<SelDefSpill> sels_spill = {
  {"kNoSpillCut",           kNoSpillCut,           colorNoCut,           "NO SELECTION CUTS APPLIED."},
  {"kContainedSpill",       kContainedSpillCut,    colorContainedCut,    "Existence of a slice in this spill passing kNuEContainedFDCut."},
  {"kFlashMatchSpill",      kFlashMatchSpillCut,   colorFlashMatchCut,   "Existence of a slice in this spill passing kSlcFlashMatchCut."},
  {"kRecoShowerFDSpillCut", kRecoShowerFDSpillCut, colorRecoShowerFDCut, "Existence of a slice in this spill passing kRecoShowerFDCut."}
};

//! @note: There are no crtsels_slice entries since CRT cuts are all at
//!        the spill level, i.e. SpillCuts.
std::vector<SelDefSpill> crtsels_spill = {
  {"kNoSpillCut",   kNoSpillCut,   colorNoCut,                "NO SELECTION CUTS APPLIED."},
  {"kCRTHitVetoFD", kCRTHitVetoFD, colorCRTHitVetoFDSpillCut, "All CRT hits within time window are NOT above 100 PE."}
};
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Interaction Types:

//! @note: kSlcIsRecoNu is A RECONSTRUCTION QUANTITY, !slc->is_clear_cosmic,
//!        which is determined at the ArtROOT level.
std::vector<SelDefSlice> types_slice = {
  {"NuTotal",         kSlcIsRecoNu, colorTotal,   "All slices with a reconstructed neutrino."},
  {"NuECC",           kRecoNuECC,   colorNuE,     "True charged-current electron neutrino (NuECC) in a slice with a reconstructed neutrino."},
  {"NuECCTrueSignal", kRecoNuECC,   colorNuE,     "True NuECC (additionally meeting containment and shower energy restrictions) in a slice with a reconstructed neutrino."},
  {"NuMuCC",          kRecoNuMuCC,  colorNuMu,    "True charged-current muon neutrino (NuMuCC) in a slice with a reconstructed neutrino."},
  {"NC",              kRecoNC,      colorNC,      "True neutral-current (NC) neutrino in a slice with a reconstructed neutrino."},
  {"NuCosmic",        kRecoCosmic,  colorCosmics, "True cosmic in a slice with a reconstructed neutrino [ALWAYS BACKGROUND]."}
};

std::vector<SelDefSpill> types_spill = {
  {"TotalSpill",  kNoSpillCut,     colorTotal,   "NO SELECTION CUTS APPLIED."},
  {"NuECCSpill",  kTrueNuECCSpillCut,  colorNuE,     "Not a cosmic/single-nu spill; first slice has e+/e- and is not NC."},
  {"NuMuCCSpill", kTrueNuMuCCSpillCut, colorNuMu,    "Not a cosmic-single-nu spill; first slice has mu+/mu- and it not NC."},
  {"NCSpill",     kTrueNCSpillCut,     colorNC,      "Not a cosmic/single-nu spill; first slice is NC."},
  {"CosmicSpill", kTrueCosmicSpillCut, colorCosmics, "Spill has no neutrinos."}
};
//! ////////////////////////////////////////////////////////////////////////////