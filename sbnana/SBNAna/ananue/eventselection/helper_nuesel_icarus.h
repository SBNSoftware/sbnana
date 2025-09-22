#pragma once 

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
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

#include <vector>
#include <algorithm>
#include <iostream>
#include <type_traits> //! required for std::is_same


#include "TStyle.h"

using namespace ana;


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


std::pair<int,int> slice_indices(int total, int nRanks, int rank) {
    // number of items per slice, some slices may have 1 extra
    int base = total / nRanks;
    int rem  = total % nRanks; // leftover items to distribute

    int start = rank * base + std::min(rank, rem);
    int end   = start + base + (rank < rem ? 1 : 0);

    // clamp to total
    end = std::min(end, total);

    return {start, end};
}


//! ////////////////////////////////////////////////////////////////////////////
//! These are examples of useful structs to making a bunch of Spectra.
//! @note In the main script, make_spectra_nuesel_icarus.C, we have vectors of,
//! e.g. plots, that sometimes consider slices and other times consider spills.
//! The templates below allow for treatment of both slices and spills and make 
//! the main script slightly cleaner.

//! Plot structure:
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
struct PlotDef
{
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  Var var = kCounting;
};
struct PlotDefSpill
{
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  SpillVar var = kSpillCounting;
};
struct PlotDefMultiVar
{
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  SpillMultiVar var = kCRTHitX;
};

//! Selection structure:
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
struct SelDef
{
  std::string label = "";
  Cut cut = kNoCut;
  int color = kBlack;
};
struct SelDefSpill
{
  std::string label = "";
  SpillCut cut = kNoSpillCut;
  int color = kBlack;
};


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

//! ////////////////////////////////////////////////////////////////////////////
//! Spectrum-building and Spectrum-saving functions:
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSpectra(
      const std::vector<PlotT>& plots,
      const std::vector<SelT>& sels,
      const std::vector<SelT>& types,
      ana::SpectrumLoaderBase& loader,
      int rank, int nRanks,
      const std::string& finname = "notCosmic",
      const double& cosmicPOT = -1.0)
{
  const unsigned int kNVar  = plots.size();
  const unsigned int kNSel  = sels.size();
  const unsigned int kNType = types.size();

  std::vector<std::vector<std::vector<ana::Spectrum*>>> specs( kNSel, 
    std::vector<std::vector<ana::Spectrum*>>( kNType, 
      std::vector<ana::Spectrum*>( kNVar, nullptr)));

  int total = kNVar * kNType * kNSel;
  auto [start, end] = slice_indices( total, nRanks, rank);

  int idx = 0;
  for (unsigned int iSel = 0; iSel < kNSel; ++iSel) {
    for (unsigned int iType = 0; iType < kNType; ++iType) {
      for (unsigned int iVar = 0; iVar < kNVar; ++iVar, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in slices

        if constexpr (std::is_same<SelT, SelDef>::value) { //! slice-level sels and types
          specs[iSel][iType][iVar] = new ana::Spectrum(
            plots[iVar].label + "_" + sels[iSel].label + "_" + types[iType].label,
            plots[iVar].bins,
            loader,
            plots[iVar].var,
            kNoSpillCut,
            sels[iSel].cut && types[iType].cut);
        }
        else if constexpr (std::is_same<SelT, SelDefSpill>::value){ //! spill-level sels and types
          specs[iSel][iType][iVar] = new ana::Spectrum(
            plots[iVar].label + "_" + sels[iSel].label + "_" + types[iType].label,
            plots[iVar].bins,
            loader,
            plots[iVar].var,
            sels[iSel].cut && types[iType].cut);
        }
        else {
          std::cerr << "[ERROR] Type SelT is neither SelDef nor SelDefSpill. Stopping program." << std::endl;
          std::exit(101);
        }

        if( finname == "cosmics") {
          std::cout << "[WARN] Fake POT scaling of cosmics!" << std::endl;

          const double specPOT( specs[iSel][iType][iVar]->POT());

          if( specPOT < std::numeric_limits<double>::epsilon()) {
            specs[iSel][iType][iVar]->OverridePOT( cosmicPOT);
          }
        }
      }
    }
  }

  return specs;
}


template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSpectraVeto(
      const std::vector<PlotT>& plots,
      const std::vector<SelT>& sels,
      const std::vector<SelT>& types,
      ana::SpectrumLoaderBase& loader,
      int rank, int nRanks,
      SpillCut kVetoCut = kNoSpillCut,
      const std::string& finname = "not_cosmic",
      const double& cosmicPOT = -1.0)
{
  const unsigned int kNVar  = plots.size();
  const unsigned int kNSel  = sels.size();
  const unsigned int kNType = types.size();

  std::vector<std::vector<std::vector<ana::Spectrum*>>> specs_veto( kNSel,
    std::vector<std::vector<ana::Spectrum*>>( kNType, 
      std::vector<ana::Spectrum*>( kNVar, nullptr)));

  int total = kNVar * kNType * kNSel;
  auto [start, end] = slice_indices( total, nRanks, rank);

  int idx = 0;
  for (unsigned int iSel = 0; iSel < kNSel; ++iSel) {
    for (unsigned int iType = 0; iType < kNType; ++iType) {
      for (unsigned int iVar = 0; iVar < kNVar; ++iVar, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in slices

        if constexpr (std::is_same<SelT, SelDef>::value) { //! slice-level sels and types
          specs_veto[iSel][iType][iVar] = new ana::Spectrum(
            plots[iVar].label + "_" + sels[iSel].label + "_" + types[iType].label + "_veto",
            plots[iVar].bins,
            loader,
            plots[iVar].var,
            kVetoCut,
            sels[iSel].cut && types[iType].cut);
        }
        else if constexpr (std::is_same<SelT, SelDefSpill>::value){ //! spill-level sels and types
          specs_veto[iSel][iType][iVar] = new ana::Spectrum(
            plots[iVar].label + "_" + sels[iSel].label + "_" + types[iType].label,
            plots[iVar].bins,
            loader,
            plots[iVar].var,
            kVetoCut && sels[iSel].cut && types[iType].cut);
        }
        else {
          std::cerr << "[ERROR] Type SelT is neither SelDef nor SelDefSpill. Stopping program." << std::endl;
          std::exit(101);
        }

        if( finname == "cosmics") {
          std::cout << "[WARN] Fake POT scaling of cosmics!" << std::endl;

          const double specPOT( specs_veto[iSel][iType][iVar]->POT());

          if( specPOT < std::numeric_limits<double>::epsilon()) {
            specs_veto[iSel][iType][iVar]->OverridePOT( cosmicPOT);
          }
        }
      }
    }
  }

  return specs_veto;
}

//! @note CRT Spectra are always for Spill level. No templating needed here.
inline std::vector<std::vector<ana::Spectrum*>> buildCRTSpectra(  
  const std::vector<PlotDefMultiVar>& plots,
  const std::vector<SelDefSpill>& sels,
  ana::SpectrumLoaderBase& loader,
  int rank, int nRanks,
  SpillCut kCRTCut = kNoSpillCut)
{
  const unsigned int kNVarCRTSpill = plots.size();
  const unsigned int kNSelCRTSpill = sels.size();

  std::vector<std::vector<ana::Spectrum*>> specs_crt( kNSelCRTSpill,
    std::vector<ana::Spectrum*>( kNVarCRTSpill, nullptr));

  int total = kNVarCRTSpill * kNSelCRTSpill;
  auto [start, end] = slice_indices( total, nRanks, rank);

  int idx = 0;
  for( unsigned int iSel = 0; iSel < kNSelCRTSpill; ++iSel) {
    for( unsigned int iVar = 0; iVar < kNVarCRTSpill; ++iVar, ++idx) {
      if( idx < start || idx >= end) continue; //! only process in slices

      specs_crt[iSel][iVar] = new ana::Spectrum(
        plots[iVar].label + "_" + sels[iSel].label + "_CRT",
        plots[iVar].bins,
        loader,
        plots[iVar].var,
        kCRTCut && sels[iSel].cut);
    }
  }

  return specs_crt;
}

//! @brief Iterate over a 3-dimensional vector of selection, interaction type,
//!        and variable-to-plot indices and save the constituent Spectra to an
//!        output ROOT file.
template <typename PlotT, typename SelT>
inline void saveSpectra( 
      TFile& fout,
      const std::vector<std::vector<std::vector<ana::Spectrum*>>>& specs,
      const std::vector<PlotT>& plots,
      const std::vector<SelT>& types,
      const std::vector<SelT>& sels,
      const int rank,
      const int nRanks,
      const std::string& vetoFlag = "") {
  const unsigned int kNVar  = plots.size();
  const unsigned int kNSel  = sels.size();
  const unsigned int kNType = types.size();

  int total = kNVar * kNType * kNSel;
  auto [start, end] = slice_indices( total, nRanks, rank);

  int idx = 0;
  for( unsigned int iSel = 0; iSel < sels.size(); ++iSel ) {
    for( unsigned int iType = 0; iType < types.size(); ++iType) {
      for( unsigned int iVar = 0; iVar < plots.size(); ++iVar, ++idx) {  
        if( idx < start || idx >= end) continue; //! only process in slices 
        ana::Spectrum *spec = specs.at(iSel).at(iType).at(iVar);

        std::string label = plots[iVar].label + "_" + sels[iSel].label + "_" + types[iType].label + "_" + vetoFlag;

        std::cout << "[INFO] Saving spectra: " << label << std::endl;
        spec->SaveTo( fout.mkdir( label.c_str()));
      }
    }
  }
}

//! @brief Iterate over a 2-dimensional vector of selection and variable-to-plot
//!        indices and save the constituent CRT Spectra to an output ROOT file.
//! @note No templating is done since CRT Spectra are always a spill-level plot.
inline void saveCRTSpectra( 
      TFile& fout,
      const std::vector<std::vector<ana::Spectrum*>>& specs_crt,
      const std::vector<PlotDefMultiVar>& plots,
      const std::vector<SelDefSpill>& sels,
      const int rank,
      const int nRanks)
{
  const unsigned int kNVar  = plots.size();
  const unsigned int kNSel  = sels.size();

  int total = kNVar * kNSel;
  auto [start, end] = slice_indices( total, nRanks, rank);

  int idx = 0;
  for (size_t iSel = 0; iSel < sels.size(); ++iSel) {
    for (size_t iVar = 0; iVar < plots.size(); ++iVar, ++idx) {
      if( idx < start || idx >= end) continue; //! only process in slices
      ana::Spectrum *spec = specs_crt.at(iSel).at(iVar);

      std::string label = plots[iVar].label + "_" + sels[iSel].label + "_CRT"; 

      std::cout << "[INFO] Saving spectra: " << label << std::endl;
      spec->SaveTo( fout.mkdir( label.c_str()));
    }
  }
}
//! ////////////////////////////////////////////////////////////////////////////