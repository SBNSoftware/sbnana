//! ////////////////////////////////////////////////////////////////////////////
//! @file: generateSpectra_icarus.h                                           //
//! @authors: Jacob Smith (smithja)                                           //
//! @email: jacob.a.smith@stonybrook.edu                                      //
//! Last edited: November 24th, 2025                                          //
//!                                                                           //
//! @details: Header file for spectrum building and saving functions          //
//! used in the ICARUS electron neutrino selection analysis.                  //
//! ////////////////////////////////////////////////////////////////////////////
#pragma once

#include <vector>      //! to use std::vector objects
#include <string>      //! to use std::string objects
#include <iostream>    //! to use std::cout
#include <type_traits> //! to use std::is_same

#include "TFile.h"       //! ROOT file I/O
#include "TTreeReader.h" //! ROOT tree reading
#include "TTree.h"       //! ROOT tree structure

//! CAFAna objects used in this file:
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

//! Specifics of the electron neutrino selection are located in this file:
#include "nuESelectionHelper_icarus.h"

#include "nuESelection_Cuts.h"

//! Imports user-defined general tools developed over a research career:
#include "auxTools.h"

using namespace ana;

//! ============================================================================
//! @section Function Declarations
//! ============================================================================

/** @fn buildSpectra()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSpectra(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT);

/** @fn buildSequentialSpectra()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types with sequential cut application.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSequentialSpectra(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT);

/** @fn buildSpectraVeto()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types. This function is virtually the same as buildSpectra(), but
 * it makes use of the CRT Veto.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSpectraVeto(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT);

/** @fn buildSequentialSpectraVeto()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types with sequential cut application. This function is virtually
 * the same as buildSequentialSpectra(), but it makes use of the CRT Veto.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSequentialSpectraVeto(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT);

/** @fn buildCRTSpectra()
 * @brief Create Spectrum objects from arrays of CRT variables and selection cuts.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 *
 * @return 2D vector of pointers for Spectra created with different variables
 * (vars) and cuts (sels)
 *
 * @note CRT Spectra are always for Spill level. No templating needed here.
*/
inline std::vector<std::vector<ana::Spectrum*>> buildCRTSpectra(
    const std::vector<PlotDefMultiVar>& vars,
    const std::vector<SelDefSpill>& sels,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks);

/** @fn saveCRTSpectra()
 * @brief Save CRT Spectrum objects to an output ROOT file.
 *
 * @param fOut: Output Tfile to which you want to save Spectra
 * @param specs_crt: 2D vector of the Spectra to save
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
*/
inline void saveCRTSpectra(
    TFile& fOut,
    const std::vector<std::vector<ana::Spectrum*>>& specs_crt,
    const std::vector<PlotDefMultiVar>& vars,
    const std::vector<SelDefSpill>& sels,
    const int rank,
    const int nRanks);

//! ============================================================================
//! @section Template and Inline Function Implementations
//! ============================================================================

/** @fn buildSpectra()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSpectra(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT) {
  const size_t kNVar  = vars.size();
  const size_t kNSel  = sels.size();
  const size_t kNType = types.size();

  std::vector<std::vector<std::vector<ana::Spectrum*>>> specs( kNVar,
    std::vector<std::vector<ana::Spectrum*>>( kNSel,
      std::vector<ana::Spectrum*>( kNType, nullptr)));

  int total = kNVar * kNSel * kNType;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for (size_t iVar = 0; iVar < kNVar; ++iVar) {
    for (size_t iSel = 0; iSel < kNSel; ++iSel) {
      for (size_t iType = 0; iType < kNType; ++iType, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in batches

        if constexpr (std::is_same<SelT, SelDefSlice>::value) { //! slice-level sels and types
          specs[iVar][iSel][iType] = new ana::Spectrum(
            "var=" + vars[iVar].label \
                  + "_unit=" + vars[iVar].unit \
                  + "_sel=" + sels[iSel].label \
                  + "_type=" + types[iType].label,
            vars[iVar].bins,
            loader,
            vars[iVar].var,
            kNoSpillCut,
            sels[iSel].cut && types[iType].cut);
        }
        else if constexpr (std::is_same<SelT, SelDefSpill>::value){ //! spill-level sels and types
          specs[iVar][iSel][iType] = new ana::Spectrum(
            "var=" + vars[iVar].label \
                  + "_unit=" + vars[iVar].unit \
                  + "_sel=" + sels[iSel].label \
                  + "_type=" + types[iType].label,
            vars[iVar].bins,
            loader,
            vars[iVar].var,
            sels[iSel].cut && types[iType].cut);
        }
        else {
          std::cerr << "[ERROR] Type SelT is neither SelDefSlice nor SelDefSpill.";
          std::cerr << " Stopping program." << std::endl;
          std::exit(101);
        }

        if( isCosmic) {
          std::cout << "[WARN] Fake POT scaling of cosmics!" << std::endl;

          const double specPOT = specs[iVar][iSel][iType]->POT();

          if( specPOT < std::numeric_limits<double>::epsilon()) {
            specs[iVar][iSel][iType]->OverridePOT( cosmicPOT);
          }
        }
      }
    }
  }

  return specs;
}


/** @fn buildSequentialSpectra()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSequentialSpectra(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT) {
  if constexpr ( not std::is_same<SelT, SelDefSlice>::value) {
    std::cerr << "[ERROR] Selection Type (SelT) is not SelDefSlice.";
    std::cerr << " Only do sequential cuts at slice level.";
    std::cerr << " Stopping program." << std::endl;
    std::exit(101);
  }

  const size_t kNVar  = vars.size();
  const size_t kNSel  = sels.size();
  const size_t kNType = types.size();

  std::vector<std::vector<std::vector<ana::Spectrum*>>> seq_specs( kNVar,
    std::vector<std::vector<ana::Spectrum*>>( kNSel,
      std::vector<ana::Spectrum*>( kNType, nullptr)));

  int total = kNVar * kNSel * kNType;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for (size_t iVar = 0; iVar < kNVar; ++iVar) {
    //! @attention: we start at index 2 in the selection loop below
    //! since the first index of the individual_cuts_slice vector will
    //! be a singular cut, e.g. kNuEContainedFDCut, and the zero-th index will
    //! always be the kNoCut case. This avoides the following error:
    //!
    //! Error in <TFile::mkdir>: An object with name <name> exists already
    //!
    //! Moreover, we initialize the first kSequentialCut to be the 1 (i.e.
    //! second) index of sels to avoid this duplication logic.
    Cut kSequentialCut = sels[1].cut; //! <-- first non-kNoCut Cut

    for (size_t iSel = 2; iSel < kNSel - 2; ++iSel) {
      kSequentialCut = kSequentialCut && sels[iSel].cut; //! add next cut

      for (size_t iType = 0; iType < kNType; ++iType, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in batches

        //! @note: sequential labeling done via iSel index; meaningful
        //!        name too long to place in Spectrum label so you must infer
        //!        what the Nth sequential cut is
        seq_specs[iVar][iSel][iType] = new ana::Spectrum(
          "var=" + vars[iVar].label \
                + "_unit=" + vars[iVar].unit \
                + "_sel=SEQ" + std::to_string( iSel) \
                + "_type=" + types[iType].label,
          vars[iVar].bins,
          loader,
          vars[iVar].var,
          kNoSpillCut,
          kSequentialCut && types[iType].cut);

        if( isCosmic) {
          std::cout << "[WARN] Fake POT scaling of cosmics!" << std::endl;

          const double seqSpecPOT = seq_specs[iVar][iSel][iType]->POT();

          if( seqSpecPOT < std::numeric_limits<double>::epsilon()) {
            seq_specs[iVar][iSel][iType]->OverridePOT( cosmicPOT);
          }
        }
      }
    }
  }

  return seq_specs;
}

/** @fn buildSpectraVeto()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types. This function is virtually the same as buildSpectra(), but
 * it makes use of the CRT Veto.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSpectraVeto(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT) {
  const size_t kNVar  = vars.size();
  const size_t kNSel  = sels.size();
  const size_t kNType = types.size();

  std::vector<std::vector<std::vector<ana::Spectrum*>>> specs_veto( kNVar,
    std::vector<std::vector<ana::Spectrum*>>( kNSel,
      std::vector<ana::Spectrum*>( kNType, nullptr)));

  int total = kNVar * kNSel * kNType;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for (size_t iVar = 0; iVar < kNVar; ++iVar) {
    for (size_t iSel = 0; iSel < kNSel; ++iSel) {
      for (size_t iType = 0; iType < kNType; ++iType, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in batches

        if constexpr (std::is_same<SelT, SelDefSlice>::value) { //! slice-level sels and types
          specs_veto[iVar][iSel][iType] = new ana::Spectrum(
            "var=" + vars[iVar].label \
                  + "_unit=" + vars[iVar].unit \
                  + "_sel=" + sels[iSel].label \
                  + "_type=" + types[iType].label \
                  + "_kCRTHitVetoFD",
            vars[iVar].bins,
            loader,
            vars[iVar].var,
            kCRTHitVetoFD,
            sels[iSel].cut && types[iType].cut);
        }
        else if constexpr (std::is_same<SelT, SelDefSpill>::value){ //! spill-level sels and types
          specs_veto[iVar][iSel][iType] = new ana::Spectrum(
            "var=" + vars[iVar].label \
                  + "_unit=" + vars[iVar].unit \
                  + "_sel=" + sels[iSel].label \
                  + "_type=" + types[iType].label \
                  + "_kCRTHitVetoFD",
            vars[iVar].bins,
            loader,
            vars[iVar].var,
            kCRTHitVetoFD && sels[iSel].cut && types[iType].cut);
        }
        else {
          std::cerr << "[ERROR] Type SelT is neither SelDefSlice nor SelDefSpill.";
          std::cerr << " Stopping program." << std::endl;
          std::exit(101);
        }

        if( isCosmic) {
          std::cout << "[WARN] Fake POT scaling of cosmics!" << std::endl;

          const double specPOT = specs_veto[iVar][iSel][iType]->POT();

          if( specPOT < std::numeric_limits<double>::epsilon()) {
            specs_veto[iVar][iSel][iType]->OverridePOT( cosmicPOT);
          }
        }
      }
    }
  }

  return specs_veto;
}

/** @fn buildSequentialSpectraVeto()
 * @brief Create Spectrum objects from arrays of variables, selection cuts, and
 * interaction types. This function is virtually the same as buildSpectra(), but
 * it makes use of the CRT Veto.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param isCosmic: Indicates if data is cosmic; default is non-cosmic
 * @param cosmicPOT: Protons on Target for cosmic data
 *
 * @return 3D vector of pointers for Spectra created with different variables
 * (vars), cuts (sels), and interaction types (types)
*/
template <typename PlotT, typename SelT>
std::vector<std::vector<std::vector<ana::Spectrum*>>> buildSequentialSpectraVeto(
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks,
    const bool isCosmic,
    const double& cosmicPOT) {
  if constexpr ( not std::is_same<SelT, SelDefSlice>::value) {
    std::cerr << "[ERROR] Selection Type (SelT) is not SelDefSlice.";
    std::cerr << " Only do sequential cuts at slice level.";
    std::cerr << " Stopping program." << std::endl;
    std::exit(101);
  }

  const size_t kNVar  = vars.size();
  const size_t kNSel  = sels.size();
  const size_t kNType = types.size();

  std::vector<std::vector<std::vector<ana::Spectrum*>>> seq_specs_veto( kNVar,
    std::vector<std::vector<ana::Spectrum*>>( kNSel,
      std::vector<ana::Spectrum*>( kNType, nullptr)));

  int total = kNVar * kNSel * kNType;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for (size_t iVar = 0; iVar < kNVar; ++iVar) {
    //! @attention: we start at index 2 in the selection loop below
    //! since the first index of the individual_cuts_slice vector will
    //! be a singular cut, e.g. kNuEContainedFDCut, and the zero-th index will
    //! always be the kNoCut case. This avoides the following error:
    //!
    //! Error in <TFile::mkdir>: An object with name <name> exists already
    //!
    //! Moreover, we initialize the first kSequentialCut to be the first
    //! index of sels to avoid this duplication logic.
    Cut kSequentialCut = sels[1].cut; //! <-- first non-kNoCut Cut

    for (size_t iSel = 2; iSel < kNSel - 2; ++iSel) {
      kSequentialCut = kSequentialCut && sels[iSel].cut; //! add next cut

      for (size_t iType = 0; iType < kNType; ++iType, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in batches

        seq_specs_veto[iVar][iSel][iType] = new ana::Spectrum(
          "var=" + vars[iVar].label \
                + "_unit=" + vars[iVar].unit \
                + "_sel=SEQ" + std::to_string( iSel) \
                + "_type=" + types[iType].label \
                + "_kCRTHitVetoFD",
          vars[iVar].bins,
          loader,
          vars[iVar].var,
          kCRTHitVetoFD,
          kSequentialCut && types[iType].cut);

        if( isCosmic) {
          std::cout << "[WARN] Fake POT scaling of cosmics!" << std::endl;

          const double seqSpecPOT = seq_specs_veto[iVar][iSel][iType]->POT();

          if( seqSpecPOT < std::numeric_limits<double>::epsilon()) {
            seq_specs_veto[iVar][iSel][iType]->OverridePOT( cosmicPOT);
          }
        }
      }
    }
  }

  return seq_specs_veto;
}

/** @fn buildCRTSpectra()
 * @brief Create Spectrum objects from arrays of CRT variables and selection cuts.
 *
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param loader: Object that actually does the filling of histograms/Spectra
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 *
 * @return 2D vector of pointers for Spectra created with different variables
 * (vars) and cuts (sels)
 *
 * @note CRT Spectra are always for Spill level. No templating needed here.
*/
inline std::vector<std::vector<ana::Spectrum*>> buildCRTSpectra(
    const std::vector<PlotDefMultiVar>& vars,
    const std::vector<SelDefSpill>& sels,
    ana::SpectrumLoaderBase& loader,
    int rank,
    int nRanks)
{
  const size_t kNVarCRTSpill = vars.size();
  const size_t kNSelCRTSpill = sels.size();

  std::vector<std::vector<ana::Spectrum*>> specs_crt( kNVarCRTSpill,
    std::vector<ana::Spectrum*>( kNSelCRTSpill, nullptr));

  int total = kNVarCRTSpill * kNSelCRTSpill;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for( size_t iVar = 0; iVar < kNVarCRTSpill; ++iVar) {
    for( size_t iSel = 0; iSel < kNSelCRTSpill; ++iSel, ++idx) {
      if( idx < start || idx >= end) continue; //! only process in batches

      specs_crt[iVar][iSel] = new ana::Spectrum(
        "var=" + vars[iVar].label \
              + "_unit=" + vars[iVar].unit \
              + "_sel=" + sels[iSel].label \
              + "_type=CRT",
        vars[iVar].bins,
        loader,
        vars[iVar].var,
        kCRTHitVetoFD && sels[iSel].cut);
    }
  }

  return specs_crt;
}

/** @fn saveSpectra()
 * @brief Save Spectrum objects to an output ROOT file.
 *
 * @param fOut: Output Tfile to which you want to save Spectra
 * @param specs: 3D vector of the Spectra to save
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param vetoFlag: Indicates if Spectra saved use CRT Veto
 * @param label: Optional custom label to prepend to directory name
*/
template <typename PlotT, typename SelT>
inline void saveSpectra(
    TFile& fOut,
    const std::vector<std::vector<std::vector<ana::Spectrum*>>>& specs,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    const int rank,
    const int nRanks,
    const std::string& vetoFlag,
    const std::string& label = "") {
  const size_t kNVar  = vars.size();
  const size_t kNSel  = sels.size();
  const size_t kNType = types.size();

  int total = kNVar * kNSel * kNType;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for( size_t iVar = 0; iVar < kNVar; ++iVar ) {
    for( size_t iSel = 0; iSel < kNSel; ++iSel) {
      for( size_t iType = 0; iType < kNType; ++iType, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in batches
        ana::Spectrum *spec = specs.at(iVar).at(iSel).at(iType);

        std::string dirLabel = "var=" + vars[iVar].label \
        + "_unit=" + vars[iVar].unit \
        + "_sel=" + sels[iSel].label \
        + "_type=" + types[iType].label\
        + "_" + vetoFlag;

        if (!label.empty()) {
          //! Prepend custom label if provided
          dirLabel = label + "_" + dirLabel;
        }

        std::cout << "[INFO] Saving spectra: " << dirLabel << std::endl;
        spec->SaveTo( fOut.mkdir( dirLabel.c_str()));
      }
    }
  }
}

/** @fn saveSequentialSpectra()
 * @brief Save Spectrum objects to an output ROOT file.
 *
 * @param fOut: Output Tfile to which you want to save Spectra
 * @param specs: 3D vector of the Spectra to save
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param types: Interaction type (e.g. NuMuCC or NC) of Spectra variables
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
 * @param vetoFlag: Indicates if Spectra saved use CRT Veto
 * @param label: Optional custom label to prepend to directory name
*/
template <typename PlotT, typename SelT>
inline void saveSequentialSpectra(
    TFile& fOut,
    const std::vector<std::vector<std::vector<ana::Spectrum*>>>& specs,
    const std::vector<PlotT>& vars,
    const std::vector<SelT>& sels,
    const std::vector<SelT>& types,
    const int rank,
    const int nRanks,
    const std::string& vetoFlag,
    const std::string& label = "") {
  const size_t kNVar  = vars.size();
  const size_t kNSel  = sels.size();
  const size_t kNType = types.size();

  int total = kNVar * kNSel * kNType;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for( size_t iVar = 0; iVar < kNVar; ++iVar ) {
    //! @attention: we start at index 2 in the selection loop below
    //! since the first index of the individual_cuts_slice vector will
    //! be a singular cut, e.g. kNuEContainedFDCut, and the zero-th index will
    //! always be the kNoCut case. This avoides the following error:
    //!
    //! Error in <TFile::mkdir>: An object with name <name> exists already

    for( size_t iSel = 2; iSel < kNSel - 2; ++iSel) {
      for( size_t iType = 0; iType < kNType; ++iType, ++idx) {
        if( idx < start || idx >= end) continue; //! only process in batches
        ana::Spectrum *spec = specs.at(iVar).at(iSel).at(iType);

        std::string dirLabel = "var=" + vars[iVar].label \
        + "_unit=" + vars[iVar].unit \
        + "_sel=SEQ" + std::to_string( iSel) \
        + "_type=" + types[iType].label\
        + "_" + vetoFlag;

        if (!label.empty()) {
          //! Prepend custom label if provided
          dirLabel = label + "_" + dirLabel;
        }

        std::cout << "[INFO] Saving spectra: " << dirLabel << std::endl;
        spec->SaveTo( fOut.mkdir( dirLabel.c_str()));
      }
    }
  }
}

/** @fn saveCRTSpectra()
 * @brief Save CRT Spectrum objects to an output ROOT file.
 *
 * @param fOut: Output Tfile to which you want to save Spectra
 * @param specs_crt: 2D vector of the Spectra to save
 * @param vars: Variables for which you want to create Spectra
 * @param sels: Cuts you want to apply to Spectra variables
 * @param rank: Process ID used when batch-processing
 * @param nRanks: Total number of process IDs when batch-processing
*/
inline void saveCRTSpectra(
    TFile& fOut,
    const std::vector<std::vector<ana::Spectrum*>>& specs_crt,
    const std::vector<PlotDefMultiVar>& vars,
    const std::vector<SelDefSpill>& sels,
    const int rank,
    const int nRanks)
{
  const size_t kNVar  = vars.size();
  const size_t kNSel  = sels.size();

  int total = kNVar * kNSel;
  auto [start, end] = sliceIndices( total, nRanks, rank);

  int idx = 0;
  for (size_t iVar = 0; iVar < kNVar; ++iVar) {
    for (size_t iSel = 0; iSel < kNSel; ++iSel, ++idx) {
      if( idx < start || idx >= end) continue; //! only process in batches
      ana::Spectrum *spec = specs_crt.at(iVar).at(iSel);

      std::string dirLabel = "var=" + vars[iVar].label \
        + "_unit=" + vars[iVar].unit \
        + "_sel=" + sels[iSel].label \
        + "_type=CRT";

      std::cout << "[INFO] Saving spectra: " << dirLabel << std::endl;
      spec->SaveTo( fOut.mkdir( dirLabel.c_str()));
    }
  }
}

/** @fn ComputeCosmicPOT()
 * @brief Computes the cosmicPOT from a given input file.
 *
 * @param inFile The input file from which to compute cosmicPOT.
 *
 * @return The computed cosmicPOT value, or -1.0 if an error occurs.
 */
double ComputeCosmicPOT(const std::string &inFile) {
  TFile f(inFile.c_str());
  if (!f.IsOpen()) {
    std::cerr << "[WARN] Could not open " << inFile << ", returning -1 for cosmicPOT" << std::endl;
    return -1.0;
  }

  TTree* myTree = (TTree*) f.Get("recTree");
  if (!myTree) {
    std::cerr << "[WARN] Could not access recTree in " << inFile << ", returning -1 for cosmicPOT" << std::endl;
    return -1.0;
  }

  //! Use these if using regular CAFs.
//  TTreeReaderValue< unsigned int> run(     cafReader, "hdr.run");
//  TTreeReaderValue< unsigned int> subrun(  cafReader, "hdr.subrun");
//  TTreeReaderValue< unsigned int> ngenevt( cafReader, "hdr.ngenevt");

  //! Use these instead if using flat CAFs.
  //! Otherwise cosmics will always have 'POT' set to zero.
  TTreeReader cafReader("recTree", &f);
  TTreeReaderValue< unsigned int> run( cafReader, "rec.hdr.run");
  TTreeReaderValue< unsigned int> subrun( cafReader, "rec.hdr.subrun");
  TTreeReaderValue< unsigned int> ngenevt( cafReader, "rec.hdr.ngenevt");

  unsigned int totalGenEvt = 0;

  //! Tally up the generated events values from each run-subrun pair that hasn't
  //! already been tallied up.
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> runSubRunMap;
  while (cafReader.Next()) {
    std::pair<unsigned int, unsigned int> p(*run, *subrun);
    if (runSubRunMap[p] == 0) { //! true if pair not logged
      ++runSubRunMap[p]; //! avoid double counting
      totalGenEvt += *ngenevt;
    }
  }

  //! Number of triggers matches the number of entries in the recTree TTree.
  //! Whereas the number of triggers from generated events is given by
  //! totalGenEvt, which we calculate above.
  const unsigned int numTriggers = myTree->GetEntries();
  const float nuPerSpill = 1.f / 21.8f; //! @note perhaps an average of nu per spill
  const float POTperSpill = 5e12f; //! @note perhaps an average POT per spill

  const float cosmicPOT = totalGenEvt * POTperSpill / (1.f - nuPerSpill);
  return cosmicPOT;
}