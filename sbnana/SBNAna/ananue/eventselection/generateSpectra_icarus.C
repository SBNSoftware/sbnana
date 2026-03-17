//! ////////////////////////////////////////////////////////////////////////////
//! @file: generateSpectra_icarus.C                                           //
//! @authors: Jacob Smith (smithja)                                           //
//! @email: jacob.a.smith@stonybrook.edu                                      //
//! Last edited: November 24th, 2025                                          //
//!                                                                           //
//! @details: Driver macro for generating Spectra in the ICARUS electron      //
//! neutrino selection. Template implementations are in the header file.      //
//! Takes a CAF and outputs a ROOT file containing the resulting histograms.  //
//! ////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TFile.h"

//! Include the header with all template implementations
#include "generateSpectra_icarus.h"

using namespace ana;

/**
 * @fn generateSpectra_icarus()
 *
 * @brief Main driver function to generate Spectra for the ICARUS electron
 * neutrino selection.
 *
 * @param inFile: Input file name (CAF/flatCAF)
 * @param outDir [optional]: Output directory for saving Spectra
 * @param input [optional] Type of input if using individual file (e.g. cosmic)
 * @param isCosmic [optional]: Indicates if data is cosmic; default is non-cosmic
 * @param spillMode [optional] Indicates processing of spills or slices
 * @param rank [optional] Rank of process (used in batch processing)
 * @param nRanks [optional] Total number of processes (used in batch processing)
*/
void generateSpectra_icarus( const std::string inFile,
                                const std::string outDir = "./outputs_spectra",
                                const std::string input = "nus",
                                const bool isCosmic = false,
                                const bool spillMode = false,
                                const int rank = 0,
                                const int nRanks = 1) {
  //! Testing sample file:
  // const std::string inFile = "refactored_g4_icarus_step2_g4step2_detsim_"
  //                               "stage0_stage1_13629324_0.flat.caf"
  //                               "-e8a8c7bb-7ad3-46b2-a5a6-2496ffc7199f.root";

  //! @note Nominal samweb definition for this electron neutrino selection:
  //! icaruspro_production_v09_89_01_02p01_2024A_ICARUS_MC_Sys_NuCos_2024A_MC_Sys_NuCos_respunCV_2ndV_flatcaf

  const std::string fOutName = outDir + "/" + input + "_" \
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

  //! Select all slices or a single slice per spill for the regular and veto
  //! Spectra.
  std::vector<std::vector<std::vector<Spectrum*>>> specs;
  std::vector<std::vector<std::vector<Spectrum*>>> specs_veto;
  std::vector<std::vector<std::vector<Spectrum*>>> seq_specs;
  std::vector<std::vector<std::vector<Spectrum*>>> seq_specs_veto;
  std::vector<std::vector<Spectrum*>> specs_crt;  //! @note CRT Spectra
                                                  //! always built/saved


  if( spillMode){ //! @attention: no sequential cuts in Spill Mode
    specs      = buildSpectra( vars_spill, sels_spill, types_spill,
                                  loader, rank, nRanks,
                                  isCosmic, cosmicPOT);

    //! Veto variations of plots use the same selection cuts as non-veto
    //! variations with the exception of the additional kCRTHitVetoFD, which
    //! is a spill-level cut, i.e. of type SpillCut. Cuts and SpillCuts are
    //! made at different levels, i.e. as different arguments in Spectrum
    //! constructors. Hence, we have different veto and non-veto variations
    //! of relevant plots.
    specs_veto = buildSpectraVeto( vars_spill, sels_spill, types_spill,
                                    loader, rank, nRanks,
                                    isCosmic, cosmicPOT);
  }
  else{ //! !spillMode means we're using slices
    specs      = buildSpectra( vars_slice, sels_slice, types_slice,
                                    loader, rank, nRanks,
                                    isCosmic, cosmicPOT);
    specs_veto = buildSpectraVeto( vars_slice, sels_slice, types_slice,
                                        loader, rank, nRanks,
                                        isCosmic, cosmicPOT);

    //! Slice-level plots can (read should) have sequentially applied cut
    //! variations of relevant plots.
    seq_specs      = buildSequentialSpectra( vars_slice, individual_cuts_slice,
                                types_slice,
                                loader, rank, nRanks,
                                isCosmic, cosmicPOT);
    seq_specs_veto = buildSequentialSpectraVeto( vars_slice,
                                individual_cuts_slice, types_slice,
                                loader, rank, nRanks,
                                isCosmic, cosmicPOT);
  }

  //! CRT plots are Spill level and are always generated regardless of
  //! if we're in Spill Mode or Slice Mode.
  specs_crt  = buildCRTSpectra( crtvars_spill, crtsels_spill,
                                  loader, rank, nRanks);

  //! This line is what actually fills all of the Spectra.
  loader.Go();

  if( spillMode) {
    saveSpectra( *fOut, specs, vars_spill, sels_spill, types_spill,
                    rank, nRanks, "noVetoApplied");

    //! No seperate function for saving Spectrum objects that were constructed
    //! differently. We just note the difference in the veto flag argument.
    saveSpectra( *fOut, specs_veto, vars_spill, sels_spill, types_spill,
                    rank, nRanks, "kCRTHitVetoFD");
  }
  else { //! !spillMode --> "slice mode"
    saveSpectra( *fOut, specs, vars_slice, sels_slice, types_slice,
                    rank, nRanks, "noVetoApplied");
    saveSpectra( *fOut, specs_veto, vars_slice, sels_slice, types_slice,
                    rank, nRanks, "kCRTHitVetoFD");

    //! Only slice mode plots should have sequential cut variations.
    saveSequentialSpectra( *fOut, seq_specs, vars_slice, individual_cuts_slice,
                              types_slice, rank, nRanks, "noVetoApplied");
    saveSequentialSpectra( *fOut, seq_specs_veto, vars_slice, individual_cuts_slice,
                              types_slice, rank, nRanks, "kCRTHitVetoFD");
  }
  saveCRTSpectra( *fOut, specs_crt, crtvars_spill, crtsels_spill, rank, nRanks);

  fOut->Write();
  fOut->Close();
}
