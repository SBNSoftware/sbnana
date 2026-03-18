//! ////////////////////////////////////////////////////////////////////////////
//! @file: runCutScan_icarus.C                                                //
//! @authors: Jacob Smith (smithja)                                           //
//! @email: jacob.a.smith@stonybrook.edu                                      //
//! Last edited: November 24th, 2025                                          //
//!                                                                           //
//! @details: Scan driver to generate Spectra for one cut-value combination  //
//! in the ICARUS electron neutrino selection cut optimization.               //
//!                                                                           //
//! Usage: cafe -bq runCutScan_icarus.C input.flat.caf ./out "ShowerEnergy" 0.30 0 1
//! ////////////////////////////////////////////////////////////////////////////
#pragma once

#include "TROOT.h"
#include "TFile.h"

#include <vector>
#include <string>
#include <iostream>
#include <limits>

//! CAFAna objects used in this file:
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

//! Imports user-defined general tools developed over a research career:
#include "auxTools.h"

//! Specifics of the electron neutrino selection are located in this file:
#include "nuESelectionHelper_icarus.h"

//! Cut definitions:
#include "nuESelection_Cuts.h"

//! Spectrum building/saving functions:
#include "generateSpectra_icarus.h"

using namespace ana;

/**
 * @fn runCutScan_icarus()
 *
 * @brief Main driver function to generate Spectra for one scan point.
 *
 * @param inFile: Input file name (CAF/flatCAF)
 * @param outDir [optional]: Output directory for saving Spectra
 * @param cutName [optional]: Name of the cut being scanned ("dEdx", etc.)
 * @param cutValue [optional]: Value for the cut being scanned
 * @param rank [optional]: Rank of process (used in batch processing)
 * @param nRanks [optional]: Total number of processes (used in batch processing)
 */
void runCutScan_icarus(const std::string &inFile,
                       const std::string &outDir = "./outputs_scan",
                       const std::string &input = "nus",
                       const bool isCosmic = false,
                       const std::string &cutName = "ShowerEnergy",
                       const double cutValue = 0.3,
                       const int rank = 0,
                       const int nRanks = 1)
{
  Color_t colorTrackScoreCut = kPink + 9; //! kCUSTOMTrackScoreCut; may be placed
                                          //! in nuESelectionHelper_icarus.h if
                                          //! it's decided to keep a cut on raw
                                          //! track score.

  std::cout << "[INFO] runCutScan_icarus starting: scanning " << cutName
            << " = " << cutValue << std::endl;

  //! Set all cuts to nominal values
  double shwEnergy_lb = NominalCuts::showerEnergy_lb;
  double showerdEdx_ub = NominalCuts::showerdEdx_ub;
  double showerConvGap_ub = NominalCuts::showerConvGap_ub;
  double showerDensity_lb = NominalCuts::showerDensity_lb;
  double nuETrackLen_ub = NominalCuts::nuETrackLen_ub;
  double flashMatch_lb = NominalCuts::flashMatch_lb;
  double flashMatch_ub = NominalCuts::flashMatch_ub;
  double barycenter_lb = NominalCuts::barycenter_lb;
  double barycenter_ub = NominalCuts::barycenter_ub;
  double trackScore_ub = NominalCuts::trackScore_ub;

  //! Override the one being scanned
  if (cutName == "ShowerEnergy") shwEnergy_lb = cutValue;
  else if (cutName == "ShowerdEdx") showerdEdx_ub = cutValue;
  else if (cutName == "ShowerConvGap") showerConvGap_ub = cutValue;
  else if (cutName == "ShowerDensity") showerDensity_lb = cutValue;
  else if (cutName == "NuETrackLen") nuETrackLen_ub = cutValue;
  else if (cutName == "FlashMatchLB") flashMatch_lb = cutValue;
  else if (cutName == "FlashMatchUB") flashMatch_ub = cutValue;
  else if (cutName == "BarycenterLB") barycenter_lb = cutValue;
  else if (cutName == "BarycenterUB") barycenter_ub = cutValue;
  else if (cutName == "TrackScore") trackScore_ub = cutValue;
  else {
    std::cout << "[ERROR] Unknown cut name: " << cutName << std::endl;
    return;
  }

  //! Build CUSTOM Cuts from the scan values
  Cut kCUSTOMShowerEnergyCut  = makeCUSTOMShowerEnergyCut(shwEnergy_lb);
  Cut kCUSTOMShowerdEdxCut    = makeCUSTOMShowerdEdxCut(showerdEdx_ub);
  Cut kCUSTOMShowerConvGapCut = makeCUSTOMShowerConvGapCut(showerConvGap_ub);
  Cut kCUSTOMShowerDensityCut = makeCUSTOMShowerDensityCut(showerDensity_lb);
  Cut kCUSTOMNuETrackLenCut   = makeCUSTOMNuETrackLenCut(nuETrackLen_ub);
  Cut kCUSTOMFlashMatchCut    = makeCUSTOMSlcFlashMatchCut(flashMatch_lb, flashMatch_ub);
  Cut kCUSTOMBarycenterCut    = makeCUSTOMBarycenterFMDeltaZCut(barycenter_lb, barycenter_ub);
  Cut kCUSTOMTrackScoreCut    = makeCUSTOMTrackScoreCut(trackScore_ub);

  //! 1. kRecoShowerCut is a rough quality cut ensuring non-zero best-plane
  //!    energy, conversion gap, and dE/dx.
  //! 2. kNueContainedFD ensures the shower is contained in fiducial volume.
  //! 3. kNueNumShowersCut ensures there's exactly one shower in the slice.
  //! 4. kNueTrackLenCut ensures longest track is less than 110 cm.

  //! Build combination CUSTOM Cuts (mirroring kRecoShowerFDCut and kFullCut).
  //! @note kCUSTOMRecoShowerFDCut does NOT include a custom cut on the raw
  //! track score just as kRecoShowerFDCut does not include such a cut...
  const Cut kCUSTOMRecoShowerFDCut = kRecoShowerCut && kNuENumShowersCut
                                    && kCUSTOMShowerEnergyCut && kCUSTOMShowerdEdxCut
                                    && kCUSTOMShowerConvGapCut && kCUSTOMShowerDensityCut
                                    && kCUSTOMNuETrackLenCut;

  //! @note ...however, we do include a cut on raw track score in kCUSTOMFullCut.
  const Cut kCUSTOMFullCut = kCUSTOMRecoShowerFDCut && kNuEContainedFDCut
                            && kCUSTOMFlashMatchCut && kCUSTOMBarycenterCut
                            && kCUSTOMTrackScoreCut;

  //! Helpers for N-1 cuts:
  const Cut kCUSTOMNonRecoShowerFD = kNuEContainedFDCut
                                      && kCUSTOMFlashMatchCut
                                      && kCUSTOMBarycenterCut
                                      && kCUSTOMTrackScoreCut;

  //! Prepare output file:
  const std::string fOutName = outDir + "/scan_" \
                                + cutName + Form("_%.3f_", cutValue) \
                                + extractFileName(inFile) \
                                + Form("_OUTPUT_SPECTRA_rank%d_of_%d.root", rank, nRanks);

  SpectrumLoader loader( inFile);

  TFile* fOut = TFile::Open(fOutName.c_str(), "RECREATE");
  if( !fOut->IsOpen()) {
    std::cout << "[ERROR] Output file (" << fOut->GetName() << ") not open, ";
    std::cout << "aborting." << std::endl;
  }
  else if( !fOut->IsWritable()) {
    std::cout << "[ERROR] Output file (" << fOut->GetName() << ") open ";
    std::cout << "but not writable, aborting." << std::endl;
  }
  else {
    fOut->cd();
  }

  std::cout << "[DEBUG] Setting gDirectory as output file ";
  std::cout << fOut->GetName() << std::endl;
  gDirectory = fOut;

  double cosmicPOT = ComputeCosmicPOT( inFile);

  //! Build the selection vector (CUSTOM Individual Cuts + CUSTOM N-1 cuts).
  std::vector<SelDefSlice> sels_scan = {
    {"kCUSTOMShowerEnergy",       kCUSTOMShowerEnergyCut,  colorShowerEnergyCut,       "Largest shower in this slice has an energy greater than " + std::to_string(shwEnergy_lb) + " MeV."},
    {"kCUSTOMShower_dEdx",        kCUSTOMShowerdEdxCut,    colorShowerdEdxCut,         "Largest shower in this slice has dEdx less than " + std::to_string(showerdEdx_ub) + " MeV per cm."},
    {"kCUSTOMShowerConvGap",      kCUSTOMShowerConvGapCut, colorShowerConvGapCut,      "Largest shower in this slice has a conversion gap less than " + std::to_string(showerConvGap_ub) + " cm."},
    {"kCUSTOMShowerDen",          kCUSTOMShowerDensityCut, colorShowerDensityCut,      "Largest shower in this slice has an energy density greater than " + std::to_string(showerDensity_lb) + "MeV per cm."},
    {"kCUSTOMShowerTrackLen",     kCUSTOMNuETrackLenCut,   colorNuETrackLenCut,        "Longest track in this slice is less than " + std::to_string(nuETrackLen_ub) + " cm."},
    {"kCUSTOMSlcFlashMatch",      kCUSTOMFlashMatchCut,    colorFlashMatchCut,         "Slice FlashMatch score is within (" + std::to_string(flashMatch_lb) + ", " + std::to_string(flashMatch_ub) + ")."},
    {"kCUSTOMBarycenterFMDeltaZ", kCUSTOMBarycenterCut,    colorBarycenterFMDeltaZCut, "Slice barycenterFM.deltaZ_Trigger within [" + std::to_string(barycenter_lb) + ", " + std::to_string(barycenter_ub) + ")."},
    {"kCUSTOMTrackScore",         kCUSTOMTrackScoreCut,    colorTrackScoreCut,         "Longest track in this slice has a track score within [0.5, " + std::to_string(trackScore_ub) + ")."}
  };

  sels_scan.push_back({"kCUSTOMN1ShowerEnergy",
                        kCUSTOMNonRecoShowerFD  && kRecoShowerCut && kNuENumShowersCut
                                                     && kCUSTOMShowerdEdxCut
                          && kCUSTOMShowerConvGapCut && kCUSTOMShowerDensityCut
                          && kCUSTOMNuETrackLenCut,
                        colorShowerEnergyCut, "CUSTOM N-1: no CUSTOM ShowerEnergy"});
  sels_scan.push_back({"kCUSTOMN1ShowerdEdx",
                        kCUSTOMNonRecoShowerFD  && kRecoShowerCut && kNuENumShowersCut
                          && kCUSTOMShowerEnergyCut
                          && kCUSTOMShowerConvGapCut && kCUSTOMShowerDensityCut
                          && kCUSTOMNuETrackLenCut,
                        colorShowerdEdxCut, "CUSTOM N-1: no CUSTOM dEdx"});
  sels_scan.push_back({"kCUSTOMN1ShowerConvGap",
                        kCUSTOMNonRecoShowerFD  && kRecoShowerCut && kNuENumShowersCut
                          && kCUSTOMShowerEnergyCut  && kCUSTOMShowerdEdxCut
                                                     && kCUSTOMShowerDensityCut
                          && kCUSTOMNuETrackLenCut,
                        colorShowerConvGapCut, "CUSTOM N-1: no CUSTOM ConvGap"});
  sels_scan.push_back({"kCUSTOMN1ShowerDensity",
                        kCUSTOMNonRecoShowerFD  && kRecoShowerCut && kNuENumShowersCut
                          && kCUSTOMShowerEnergyCut  && kCUSTOMShowerdEdxCut
                          && kCUSTOMShowerConvGapCut
                          && kCUSTOMNuETrackLenCut,
                        colorShowerDensityCut, "CUSTOM N-1: no CUSTOM Density"});
  sels_scan.push_back({"kCUSTOMN1TrackLen",
                        kCUSTOMNonRecoShowerFD  && kRecoShowerCut && kNuENumShowersCut
                          && kCUSTOMShowerEnergyCut  && kCUSTOMShowerdEdxCut
                          && kCUSTOMShowerConvGapCut && kCUSTOMShowerDensityCut,
                        colorNuETrackLenCut, "CUSTOM N-1: no CUSTOM TrackLen"});

  sels_scan.push_back({"kCUSTOMN1FlashMatch",
                        kCUSTOMRecoShowerFDCut && kNuEContainedFDCut
                                                  && kCUSTOMBarycenterCut
                          && kCUSTOMTrackScoreCut,
                        colorFlashMatchCut, "CUSTOM N-1: no CUSTOM FlashMatch"});
  sels_scan.push_back({"kCUSTOMN1Barycenter",
                        kCUSTOMRecoShowerFDCut && kNuEContainedFDCut
                          && kCUSTOMFlashMatchCut
                          && kCUSTOMTrackScoreCut,
                        colorBarycenterFMDeltaZCut, "CUSTOM N-1: no CUSTOM Barycenter"});
  sels_scan.push_back({"kCUSTOMN1TrackScore",
                        kCUSTOMRecoShowerFDCut && kNuEContainedFDCut
                          && kCUSTOMFlashMatchCut && kCUSTOMBarycenterCut,
                        colorTrackScoreCut, "CUSTOM N-1: no CUSTOM TrackScore"});

  //! Build spectra (no veto and with CRT veto) - following generateSpectra_icarus.C pattern
  std::vector<std::vector<std::vector<Spectrum*>>> specs =
    buildSpectra(vars_slice, sels_scan, types_slice, loader,
                  rank, nRanks, isCosmic, cosmicPOT);

  std::vector<std::vector<std::vector<Spectrum*>>> specs_veto =
    buildSpectraVeto(vars_slice, sels_scan, types_slice, loader,
                      rank, nRanks, isCosmic, cosmicPOT);

  //! This line is what actually fills all of the Spectra
  loader.Go();

  //! Save spectra with scan label in directory names
  std::string scanLabel = "scan=" + cutName + Form("_%.3f", cutValue);
  saveSpectra( *fOut, specs, vars_slice, sels_scan, types_slice,
              rank, nRanks, "noVetoApplied", scanLabel);
  saveSpectra( *fOut, specs_veto, vars_slice, sels_scan, types_slice,
              rank, nRanks, "kCRTHitVetoFD", scanLabel);

  fOut->Write();
  fOut->Close();

  std::cout << "[INFO] Finished scan job." << std::endl;
  std::cout << "[INFO] Output: " << fOutName << std::endl;
}