#pragma once 

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/NueCuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202401.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"

#include "TStyle.h"

using namespace ana;


const int limitN1 = 13; // N1 cuts start at 13 position in the sels object

//! Alias colors:
Color_t color_eff   = kGreen+2;  //! efficiency 
Color_t color_pur   = kOrange+2; //! purity
Color_t color_nue   = kBlue-7;   //! electron neutrino
Color_t color_numu  = kGreen+1;  //! muon neutrino
Color_t color_nc    = kMagenta;  //! neutral current interactions
Color_t color_cos   = kGray+2;   //! cosmics
Color_t color_other = kOrange+8;

//! Alias line styles:
Style_t line_nue    = kSolid;  //! electron neutrino
Style_t line_numu   = kSolid;  //! muon neutrino
Style_t line_nc     = kDashed; //! neutral current interactions
Style_t line_cos    = kSolid;  //! cosmics
Style_t line_other  = kDotted;



//! ////////////////////////////////////////////////////////////////////////////
//! These are examples of useful structs to making a bunch of Spectra

//! Plot structure:
struct PlotDef
{
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  Var var = kCounting;
};
struct PlotDefSpill
{
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  SpillVar var = kSpillCounting;
};
struct PlotDefMultiVar
{
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  SpillMultiVar var = kCRTHitX;
};

//! Selection structure:
struct SelDef
{
  std::string suffix = "";
  std::string label = "";
  Cut cut = kNoCut;
  int color = kBlack;
};
struct SelDefSpill
{
  std::string suffix = "";
  std::string label = "";
  SpillCut cut = kNoSpillCut;
  int color = kBlack;
};
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Define variables, cuts, and binnings:
//! @todo: give units for all relevant quantities, especially binnings

//! Variables:
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
const Var kBarycenterFM([](const caf::SRSliceProxy *slc) {
  return slc->barycenterFM.deltaZ_Trigger;
});

//! Cuts:
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
const Cut kBarycenterFMFDCut([](const caf::SRSliceProxy *slc) {
  return !std::isnan(slc->barycenterFM.deltaZ_Trigger) && 
          slc->barycenterFM.deltaZ_Trigger >= 0 && 
          slc->barycenterFM.deltaZ_Trigger < 100;
});

const Cut kTotal      = kNoCut;
const Cut kNueCC      = kIsNue && !kIsNC;  //! charged-current electron neutrino 
const Cut kNumuCC     = kIsNumu && !kIsNC; //! charged-current muon neutrino
const Cut kNC         = kIsNC;             //! neutral current interactions
const Cut kIsCosmic   = !kHasNu;           //! no neutrinos? then it's a cosmic

//! Largest shower energy > 200 MeV
const Cut kTrueFVFD = kTrueFiducialVolumeFDCryo1 || kTrueFiducialVolumeFDCryo2;
const Cut kRecoFVFD = kFiducialVolumeFDCryo1 || kFiducialVolumeFDCryo2;

//! Step by step cuts
const Cut kRecoCut   = kRecoShowerFD;
const Cut kContained = kNueContainedFD;
const Cut kFullCut   = kNueContainedFD && kNueFlashScoreFDCut && kRecoShowerFD && kBarycenterFMFDCut;

//! N-1 cuts apply all cuts except one named after 'kN1'
const Cut kN1Contained  = kNueFlashScoreFDCut && kRecoShowerFD && kBarycenterFMFDCut;
const Cut kN1Flash      = kContained && kRecoShowerFD && kBarycenterFMFDCut;
const Cut kN1Reco       = kContained && kNueFlashScoreFDCut && kBarycenterFMFDCut;
const Cut kN1RecoShower = kContained && kNueFlashScoreFDCut && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut && kBarycenterFMFDCut;
const Cut kN1NumShowers = kContained && kNueFlashScoreFDCut && kRecoShower && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut && kBarycenterFMFDCut;
const Cut kN1Dedx       = kContained && kNueFlashScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut && kBarycenterFMFDCut;
const Cut kN1ConvGap    = kContained && kNueFlashScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut && kBarycenterFMFDCut;
const Cut kN1TrkLen     = kContained && kNueFlashScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kShowerDensityCut && kShowerEnergyCut && kBarycenterFMFDCut;
const Cut kN1Density    = kContained && kNueFlashScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerEnergyCut && kBarycenterFMFDCut;
const Cut kN1Energy     = kContained && kNueFlashScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kBarycenterFMFDCut;
const Cut kN1Barycenter = kContained && kNueFlashScoreFDCut && kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut;

//! Spill Cuts:
const SpillCut kContainedSpill([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kContained(&slc))
      ++counter;
  }

  return counter;
});

const SpillCut kFlashMatchSpillCut([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kSlcFlashMatchCut(&slc))
      ++counter;
  }

  return counter;
});

const SpillCut kNuScoreSpillCut([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kSlcNuScore(&slc))
      ++counter;
  }

  return counter;
});

const SpillCut kRecoSpillCut([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kRecoCut(&slc))
      ++counter;
  }

  return counter;
});

const SpillCut kFullSpillCut([](const caf::SRSpillProxy* sr) {

  unsigned int counter(0);
  for (auto const& slc : sr->slc) {
    if (slc.tmatch.index < 0)
      continue;
    if (kFullCut(&slc))
      ++counter;
  }

  return counter;
});

const SpillCut kNueCCSpill    = kIsNueSpill && !kIsNCSpill;
const SpillCut kNumuCCSpill   = kIsNumuSpill && !kIsNCSpill;
const SpillCut kNCSpill       = kIsNCSpill;
const SpillCut kTotalSpill    = kNoSpillCut;
const SpillCut kIsCosmicSpill = kIsCosmicSpill;


//! Binnings:
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
const Binning kEnergyBinning    = Binning::Simple(40,0.,3.); // to define
const Binning kDedxBinning      = Binning::Simple(40,0.,10); // to define
const Binning kGapBinning       = Binning::Simple(40,0.,10);
const Binning kDensityBinning   = Binning::Simple(50,0.,10);
const Binning kOpenAngleBinning = Binning::Simple(60,0.,1.5);
const Binning kLengthBinning    = Binning::Simple(40,0.,200);
const Binning kPEBinning        = Binning::Simple(60,0.,600);
const Binning kTimeBinning      = Binning::Simple(155,-1550.,1550.);
const Binning kFlashBinning     = Binning::Simple(40,-6.f,34.f);
const Binning kShowerEnergyBinning  = Binning::Simple( 20, 0., 1000.); // MeV
const Binning kBarycenterFMBinning  = Binning::Simple( 41,-5., 200.);
const Binning kFMScoreBinning       = Binning::Simple( 30,  0.,   15.);
const Binning kFMTimeBinning        = Binning::Simple( 30,  -1.,   2.);
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Plots
std::vector<PlotDef> plots_slice = {
  {"count",      "Number of slices",             Binning::Simple(3,0,3), kCounting},
  {"openangle",  "Opening angle",                kOpenAngleBinning,      kRecoShower_OpenAngle},
  {"startx",     "Shower start position X (cm)", kPositionXFDBinning,    kRecoShower_StartX},
  {"starty",     "Shower start position Y (cm)", kPositionYFDBinning,    kRecoShower_StartY},
  {"startz",     "Shower start position Z (cm)", kPositionZFDBinning,    kRecoShower_StartZ},
  {"endx",       "Shower end position X (cm)",   kPositionXFDBinning,    kRecoShower_EndX},
  {"endy",       "Shower end position Y (cm)",   kPositionYFDBinning,    kRecoShower_EndY},
  {"endz",       "Shower end position Z (cm)",   kPositionZFDBinning,    kRecoShower_EndZ},
  {"vtxx",       "Slice vertex X (cm)",          kPositionXFDBinning,    kSlcVtxX},
  {"vtxy",       "Slice vertes Y (cm)",          kPositionYFDBinning,    kSlcVtxY},
  {"vtxz",       "Slice vertex Z (cm)",          kPositionZFDBinning,    kSlcVtxZ},
  {"conversion", "Conversion gap (cm)",          kGapBinning,            kRecoShower_ConversionGap},
  {"bestdedx",   "Best plane dEdx (MeV)",        kDedxBinning,           kRecoShower_BestdEdx},  
  {"bestenergy", "Best plane energy (MeV)",      kEnergyBinning,         kRecoShower_BestEnergy},  
  {"density",    "Shower density (MeV/cm)",      kDensityBinning,        kRecoShower_Density},
  {"energy",     "Shower energy (MeV)",          kShowerEnergyBinning,   kRecoShower_Energy},
  {"lengthshw",  "Shower length (cm)",           kLengthBinning,         kRecoShower_Length},
  {"lengthtrk",  "Track length (cm)",            kLengthBinning,         kLongestTrackLength},
  {"baryfm",     "Barycenter Delta Z Trigger",   kBarycenterFMBinning,   kBarycenterFM},
  {"flashscore", "Flash score",                  kFlashBinning,          kSlcFlashScore},
  {"truthenergy","True Energy (GeV)",            kLowEnergyGeVBinning,   kTruthEnergy}
};

std::vector<PlotDefSpill> plots_spill = {
  {"count",        "Number of spills",           Binning::Simple(3,0,3), kSpillCounting},
  {"nuenergy",     "True Neutrino Energy (GeV)", kLowEnergyGeVBinning,   kTruthNuEnergy},
  {"leptonenergy", "True Lepton Energy (GeV)",   kLowEnergyGeVBinning,   kTruthLeptonEnergy}
};

std::vector<PlotDefMultiVar> crtplots_spill = {
  {"crtx",    "CRT Hit Position X (cm)", kCRTXFDBinning, kCRTHitX},
  {"crty",    "CRT Hit Position Y (cm)", kCRTYFDBinning, kCRTHitY},
  {"crtz",    "CRT Hit Position Z (cm)", kCRTZFDBinning, kCRTHitZ},
  {"crttime", "CRT Hit Time (#mus)",     kTimeBinning,   kCRTHitTimeFD},
  {"crtpe",   "CRT PE",                  kPEBinning,     kCRTHitPE}
};
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Interaction types:
std::vector<SelDef> types_slice ={
  {"nuetruecont",   "NuE CC True Contained",        kSlcIsRecoNu && kNueCC && kTrueContainedFD,                                         color_nue},
  {"nuetruefv",     "NuE CC True FV",               kSlcIsRecoNu && kNueCC && kTrueFVFD,                                                color_nue},
  {"nuetruerecofv", "NuE CC True and Reco FV",      kSlcIsRecoNu && kNueCC && kTrueFVFD && kRecoFVFD,                                   color_nue},
  {"nuetruecontfv", "NuE CC True Contained and FV", kSlcIsRecoNu && kNueCC && kTrueContainedFD && kTrueFVFD,                            color_nue},
  {"nuemess",       "NuE CC Mess",                  kSlcIsRecoNu && kNueCC && kTrueContainedFD && kContained && kTrueFVFD && kRecoFVFD, color_nue},
  {"nue",           "NuE CC",                       kSlcIsRecoNu && kNueCC,                                                             color_nue},
  {"numu",          "NuMu CC",                      kSlcIsRecoNu && kNumuCC,                                                            color_numu},
  {"nunc",          "NC",                           kSlcIsRecoNu && kNC,                                                                color_nc},
  {"total",         "Total",                        kSlcIsRecoNu,                                                                       color_other},
  {"cosmic",        "Cosmic",                       kSlcIsRecoNu && kIsCosmic,                                                          color_cos}
};

//! @note kSlcIsRecoNu is already included in the spill cuts from sels_spill
std::vector<SelDefSpill> types_spill = {
  {"nue",    "NuE CC",  kNueCCSpill,      color_nue},
  {"numu",   "NuMu CC", kNumuCCSpill,     color_numu},
  {"nunc",   "NC",      kNCSpill,         color_nc},
  {"total",  "Total",   kTotalSpill,      color_other},
  {"cosmic", "Cosmic",  kIsCosmicSpill,   color_cos}
};

//! ////////////////////////////////////////////////////////////////////////////
//! Cuts:
std::vector<SelDef> sels_slice = {
  {"nocut",      "No cut",                    kNoCut,            kBlack},
  {"cont",       "Containment",               kContained,        kBlack}, //! reco fiducial volume containment
  {"flash",      "Flash score",               kSlcFlashMatchCut, kBlack},
  //! Subcuts that go into the fiducial volume reco cut:
  {"recoshw",    "Reconstructed shower",      kRecoShower,       kBlack},
  {"nshws",      "Number of showers",         kNueNumShowersCut, kBlack},
  {"dedx",       "Shower dE/dx",              kShowerdEdxCut,    kBlack},
  {"convgap",    "Conversion gap",            kShowerConvGapCut, kBlack},
  {"trklen",     "Track lenght",              kNueTrackLenCut,   kBlack},
  {"density",    "Shower density",            kShowerDensityCut, kBlack},
  {"energy",     "Shower energy",             kShowerEnergyCut,  kBlack},
  {"recocut",    "Reconstruction (all)",      kRecoCut,          kBlack},
  {"barycenter", "Barycenter",                kBarycenterFMFDCut,kBlack},
  {"everything", "Full selection",            kFullCut,          kBlack},
  //! N-1 cuts:
  {"N1cont",       "N1 Containment",          kN1Contained,  kBlack},
  {"N1flash",      "N1 Flash score",          kN1Flash,      kBlack},
  {"N1recoshw",    "N1 Reconstructed shower", kN1RecoShower, kBlack},
  {"N1nshws",      "N1 Number of showers",    kN1NumShowers, kBlack},
  {"N1dedx",       "N1 Shower dE/dx",         kN1Dedx,       kBlack},
  {"N1convgap",    "N1 Conversion gap",       kN1ConvGap,    kBlack},
  {"N1trklen",     "N1 Track length",         kN1TrkLen,     kBlack},
  {"N1density",    "N1 Shower density",       kN1Density,    kBlack},
  {"N1energy",     "N1 Shower energy",        kN1Energy,     kBlack},
  {"N1barycenter", "N1 Barycenter",           kN1Barycenter, kBlack},
  {"N1recocut",    "N1 Reconstruction (all)", kN1Reco,       kBlack}
};

std::vector<SelDefSpill> sels_spill = {
  {"nospillcut",      "No spill cut",               kNoSpillCut,         kBlack},
  {"contspill",       "Containment spill",          kContainedSpill,     kBlack},
  {"flashspill",      "Flash score spill",          kFlashMatchSpillCut, kBlack},
  {"recocutspill",    "Reconstruction (all) spill", kRecoSpillCut,       kBlack},
  {"everythingspill", "Full selection spill",       kFullSpillCut,       kBlack}
};

std::vector<SelDefSpill> crtsels_spill = {
  {"nocut_spill", "All Slices", kNoSpillCut,   kBlack},
  {"crtveto_spill", "CRTVeto",  kCRTHitVetoFD, kBlue}
};
//! ////////////////////////////////////////////////////////////////////////////