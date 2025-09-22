//! ////////////////////////////////////////////////////////////////////////////
//! Authors: Diana Patricia Mendez, Jacob Smith (smithja)                     //
//! Contact: jacob.a.smith@stonybrook.edu                                     //
//! Last edited: September 22nd, 2025                                         //   
//!                                                                           // 
//! Makes Spectra of the variables defined in helper with specific            // 
//! cuts applied to the sample(s).                                            //
//! Takes a CAF and outputs a ROOT file containing the resulting histograms.  //
//! ////////////////////////////////////////////////////////////////////////////

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "helper_nuesel_icarus.h" //! @note contains Spectrum creation templates
#include "tools_nuesel_icarus.h"

#include "TFile.h"
#include "TTreeReader.h"

#include <filesystem>
#include <iostream>
#include <string>

using namespace ana;

//! ////////////////////////////////////////////////////////////////////////////
//! Binnings:
//! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
const Binning kEnergyBinning        = Binning::Simple(40,0.,3.); //! @todo to define
const Binning kDedxBinning          = Binning::Simple(40,0.,10); //! @todo to define
const Binning kGapBinning           = Binning::Simple(40,0.,10);
const Binning kDensityBinning       = Binning::Simple(50,0.,10);
const Binning kOpenAngleBinning     = Binning::Simple(60,0.,1.5);
const Binning kLengthBinning        = Binning::Simple(40,0.,200);
const Binning kPEBinning            = Binning::Simple(60,0.,600);
const Binning kTimeBinning          = Binning::Simple(155,-1550.,1550.);
const Binning kFlashBinning         = Binning::Simple(40,-6.f,34.f);
const Binning kShowerEnergyBinning  = Binning::Simple( 20, 0., 1000.); // MeV
const Binning kBarycenterFMBinning  = Binning::Simple( 41,-5., 200.);
const Binning kFMScoreBinning       = Binning::Simple( 30,  0.,   15.);
const Binning kFMTimeBinning        = Binning::Simple( 30,  -1.,   2.);
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Plots
std::vector<PlotDef> plots_slice = {
  {"NumberOfSlices",                     Binning::Simple(3,0,3), kCounting},
  {"OpeningAngle",                       kOpenAngleBinning,      kRecoShower_OpenAngle},
  {"RecoShowerStartX_cm",                kPositionXFDBinning,    kRecoShower_StartX},
  {"RecoShowerStartY_cm",                kPositionYFDBinning,    kRecoShower_StartY},
  {"RecoShowerStartZ_cm",                kPositionZFDBinning,    kRecoShower_StartZ},
  {"RecoShowerEndX_cm",                  kPositionXFDBinning,    kRecoShower_EndX},
  {"RecoShowerEndY_cm",                  kPositionYFDBinning,    kRecoShower_EndY},
  {"RecoShowerEndZ_cm",                  kPositionZFDBinning,    kRecoShower_EndZ},
  {"SliceVertexX_cm",                    kPositionXFDBinning,    kSlcVtxX},
  {"SliceVertexY_cm",                    kPositionYFDBinning,    kSlcVtxY},
  {"SliceVertexZ_cm",                    kPositionZFDBinning,    kSlcVtxZ},
  {"RecoShowerConversionGap_cm",         kGapBinning,            kRecoShower_ConversionGap},
  {"RecoShowerBestPlane_dEdx_MeV",       kDedxBinning,           kRecoShower_BestdEdx},  
  {"RecoShowerBestPlaneEnergy_MeV",      kEnergyBinning,         kRecoShower_BestEnergy},  
  {"RecoShowerShowerDensity_MeV_per_cm", kDensityBinning,        kRecoShower_Density},
  {"RecoShowerShowerEnergy_MeV",         kShowerEnergyBinning,   kRecoShower_Energy},
  {"RecoShowerShowerLength_cm",          kLengthBinning,         kRecoShower_Length},
  {"LongestTrackLength_cm",              kLengthBinning,         kLongestTrackLength},
  {"BarycenterDeltaZTrigger",            kBarycenterFMBinning,   kBarycenterFM},
  {"FlashScore_SLICE",                   kFlashBinning,          kSlcFlashScore},
  {"TrueEnergy_GeV",                     kLowEnergyGeVBinning,   kTruthEnergy}
};

std::vector<PlotDefSpill> plots_spill = {
  {"NumberOfSpills",         Binning::Simple(3,0,3), kSpillCounting},
  {"TrueNeutrinoEnergy_GeV", kLowEnergyGeVBinning,   kTruthNuEnergy},
  {"TrueLeptonEnergy_GeV",   kLowEnergyGeVBinning,   kTruthLeptonEnergy}
};

std::vector<PlotDefMultiVar> crtplots_spill = {
  {"CRTHitX_cm",      kCRTXFDBinning, kCRTHitX},
  {"CRTHitY_cm",      kCRTYFDBinning, kCRTHitY},
  {"CRTHitZ_cm",      kCRTZFDBinning, kCRTHitZ},
  {"CRTHitTimeMus",   kTimeBinning,   kCRTHitTimeFD},
  {"CRTPE",           kPEBinning,     kCRTHitPE}
};
//! ////////////////////////////////////////////////////////////////////////////


//! ////////////////////////////////////////////////////////////////////////////
//! Interaction types:
std::vector<SelDef> types_slice = {
  {"NuECCTrueContained",      kSlcIsRecoNu && kNueCC && kTrueContainedFD,                                         color_nue},
  {"NuECCTrueFV",             kSlcIsRecoNu && kNueCC && kTrueFVFD,                                                color_nue},
  {"NuECCTrueAndRecoFV",      kSlcIsRecoNu && kNueCC && kTrueFVFD && kRecoFVFD,                                   color_nue},
  {"NuECCTrueContainedAndFV", kSlcIsRecoNu && kNueCC && kTrueContainedFD && kTrueFVFD,                            color_nue},
  {"NuECCMess",               kSlcIsRecoNu && kNueCC && kTrueContainedFD && kContained && kTrueFVFD && kRecoFVFD, color_nue},
  {"NuECC",                   kSlcIsRecoNu && kNueCC,                                                             color_nue},
  {"NuMuCC",                  kSlcIsRecoNu && kNumuCC,                                                            color_numu},
  {"NC",                      kSlcIsRecoNu && kNC,                                                                color_nc},
  {"Total",                   kSlcIsRecoNu,                                                                       color_other},
  {"Cosmic",                  kSlcIsRecoNu && kIsCosmic,                                                          color_cos}
};

//! @note kSlcIsRecoNu is already included in the spill cuts from sels_spill
std::vector<SelDefSpill> types_spill = {
  {"NuECC",  kNueCCSpill,      color_nue},
  {"NuMuCC", kNumuCCSpill,     color_numu},
  {"NC",     kNCSpill,         color_nc},
  {"Total",  kTotalSpill,      color_other},
  {"Cosmic", kIsCosmicSpill,   color_cos}
};

//! ////////////////////////////////////////////////////////////////////////////
//! Cuts:
std::vector<SelDef> sels_slice = {
  {"NoCut",       kNoCut,            kBlack},
  {"Containment", kContained,        kBlack},
  {"FlashScore",  kSlcFlashMatchCut, kBlack},
  //! Subcuts that go into the fiducial volume reco cut:
  {"RecoShower",      kRecoShower,       kBlack},
  {"NumberOfShowers", kNueNumShowersCut, kBlack},
  {"Shower_dEdx",     kShowerdEdxCut,    kBlack},
  {"ConversionGap",   kShowerConvGapCut, kBlack},
  {"TrackLength",     kNueTrackLenCut,   kBlack},
  {"ShowerDensity",   kShowerDensityCut, kBlack},
  {"ShowerEnergy",    kShowerEnergyCut,  kBlack},
  {"RecoAll",         kRecoCut,          kBlack},
  {"Barycenter",      kBarycenterFMFDCut,kBlack},
  {"FullSelection",   kFullCut,          kBlack},
  //! N-1 cuts:
  {"N-1_Containment",     kN1Contained,  kBlack},
  {"N-1_FlashScore",      kN1Flash,      kBlack},
  {"N-1_RecoShower",      kN1RecoShower, kBlack},
  {"N-1_NumberOfShowers", kN1NumShowers, kBlack},
  {"N-1_Shower_dEdx",     kN1Dedx,       kBlack},
  {"N-1_ConversionGap",   kN1ConvGap,    kBlack},
  {"N-1_TrackLength",     kN1TrkLen,     kBlack},
  {"N-1_ShowerDensity",   kN1Density,    kBlack},
  {"N-1_ShowerEnergy",    kN1Energy,     kBlack},
  {"N-1_Barycenter",      kN1Barycenter, kBlack},
  {"N-1_RecoAll",         kN1Reco,       kBlack}
};

std::vector<SelDefSpill> sels_spill = {
  {"NoSpillCut",         kNoSpillCut,         kBlack},
  {"ContainmentSpill",   kContainedSpill,     kBlack},
  {"FlashScoreSpill",    kFlashMatchSpillCut, kBlack},
  {"RecoAllSpill",       kRecoSpillCut,       kBlack},
  {"FullSelectionSpill", kFullSpillCut,       kBlack}
};

std::vector<SelDefSpill> crtsels_spill = {
  {"AllSlices", kNoSpillCut,   kBlack},
  {"CRTVeto",   kCRTHitVetoFD, kBlue}
};
//! ////////////////////////////////////////////////////////////////////////////


std::string removeFileExtension(const std::string& filename) {
    size_t lastDotPos = filename.rfind('.');
    if (lastDotPos != std::string::npos) {
        return filename.substr(0, lastDotPos);
    }
    return filename; // No extension found
}

void make_spectra_nuesel_icarus( std::string finname = "nus", 
                                 int setno = 1, 
                                 bool selspill = false,
                                 int rank = 0,
                                 int nRanks = 1) {
  std::string settag = std::to_string(setno);
  std::string foutdir = "outputs/";

  const std::string finsample = "refactored_g4_icarus_step2_g4step2_detsim_"
                                "stage0_stage1_13629324_0.flat.caf"
                                "-e8a8c7bb-7ad3-46b2-a5a6-2496ffc7199f.root"; 
//  const std::string finsample = "icaruspro_production_v09_89_01_01_2024A_ICARUS_MC_Sys_NuCos_2024A_MC_Sys_NuCos_CV_flatcaf";

  //! Create the output filename with the "spill" tag if setno is declared.
  //! @note Since setno has a default value of 1 in the make_spectra_nuesel_icarus()
  //!       function, you will likely have to delete setno to use set the output
  //!       file name to have the "slice" tag.
  const std::string foutname = foutdir + removeFileExtension(finsample) + 
      TString::Format("_OUTPUT_SPECTRA_rank%d_of_%d.root", rank, nRanks).Data();

  SpectrumLoader loader( finsample);

  TFile* fout = TFile::Open(foutname.c_str(), "RECREATE");
  if( !fout->IsOpen()) {
    std::cout << "[ERROR] Output file (" << fout->GetName() << ") not open, aborting." << std::endl;
  }
  else if( !fout->IsWritable()) {
    std::cout << "[ERROR] Output file (" << fout->GetName() << ") open but not writable, aborting." << std::endl;
  }
  else {
    fout->cd();
  }

  std::cout << "[DEBUG] Setting gDirectory as output file " << fout->GetName() << std::endl;
  gDirectory = fout;

  //////////////////////////////////////////////////////////////////////////////
  //! Pseudo-scale Cosmics:
  TFile f( finsample.c_str());
  TTree* myTree = (TTree*) f.Get("recTree");
  TTreeReader cafReader("recTree", &f);

  //! Use these if using regular CAFs.
//  TTreeReaderValue< unsigned int> run(     cafReader, "hdr.run");
//  TTreeReaderValue< unsigned int> subrun(  cafReader, "hdr.subrun");
//  TTreeReaderValue< unsigned int> ngenevt( cafReader, "hdr.ngenevt");

  //! Use these instead if using flat CAFs.
  //! Otherwise cosmics will always have 'POT' set to zero.
  TTreeReaderValue< unsigned int> run(     cafReader, "rec.hdr.run");
  TTreeReaderValue< unsigned int> subrun(  cafReader, "rec.hdr.subrun");
  TTreeReaderValue< unsigned int> ngenevt( cafReader, "rec.hdr.ngenevt");

  unsigned int totalGenEvt(0);

  //! Tally up the generated events values from each run-subrun pair that hasn't
  //! already been tallied up.
  std::map<std::pair< unsigned int, unsigned int>, unsigned int> runSubRunMap;
  while( cafReader.Next()) {
    std::pair< unsigned int, unsigned int> thisPair(*run, *subrun);
    if( runSubRunMap[thisPair] == 0) { //! evaluates to true if pair not logged
      ++runSubRunMap[thisPair]; //! avoids double-counting
      totalGenEvt += *ngenevt;
    }
  }

  //! Number of triggers matches the number of entries in the recTree TTree.
  //! Whereas the number of triggers from generated events is given by
  //! totalGenEvt, which we calculate above.
  const unsigned int numTriggers( myTree->GetEntries());


  std::cout << "Num Triggers: " << numTriggers << std::endl;
  std::cout << "Num Triggers from generated events: " << totalGenEvt << std::endl;

  const float nuPerSpill( 1.f / 21.8); //! @note perhaps some average number of
                                       //!       neutrinos per spill
  const float POTperSpill( 5e12); //! @note perhaps some average of POT per spill

  //! Weight the cosmic POT by the inverse of (1 - nuPerSpil)
  const float cosmicPOT( totalGenEvt * POTperSpill / (1.f - nuPerSpill));
  const float spills_per_cos_trig = (cosmicPOT / numTriggers) / POTperSpill; 

  std::cout << "cosmicPOT: " << cosmicPOT << std::endl;
  std::cout << "Spills per cosmic Trigger: " << spills_per_cos_trig << std::endl;
  //////////////////////////////////////////////////////////////////////////////

  //! Select all slices or a single slice per spill for the regular and veto
  //! Spectra. Most of this is abstracted away to helper_nuesel_icarus.h 
  std::vector<std::vector<std::vector<Spectrum*>>> specs;
  std::vector<std::vector<std::vector<Spectrum*>>> specs_veto;
  std::vector<std::vector<Spectrum*>> specs_crt;  //! @note CRT Spectra always built/saved


  if( selspill){
    specs      = buildSpectra( plots_spill, sels_spill, types_spill, loader, rank, nRanks);
    specs_veto = buildSpectraVeto( plots_spill, sels_spill, types_spill, loader, rank, nRanks);
  }
  else{ //! using slices
    specs      = buildSpectra( plots_slice, sels_slice, types_slice, loader, rank, nRanks);
    specs_veto = buildSpectraVeto( plots_slice, sels_slice, types_slice, loader, rank, nRanks);
  }
  specs_crt  = buildCRTSpectra( crtplots_spill, crtsels_spill, loader, rank, nRanks, kCRTHitVetoFD);

  //! This line is what actually fills all of the Spectra.
  loader.Go();

  if( selspill){
    saveSpectra( *fout, specs, plots_spill, sels_spill, types_spill, rank, nRanks, "noVetoApplied");
    saveSpectra( *fout, specs_veto, plots_spill, sels_spill, types_spill, rank, nRanks, "kCRTHitVetoFD");
  }
  else{
    saveSpectra( *fout, specs, plots_slice, sels_slice, types_slice, rank, nRanks, "noVetoApplied");
    saveSpectra( *fout, specs_veto, plots_slice, sels_slice, types_slice, rank, nRanks, "kCRTHitVetoFD");
  }
  saveCRTSpectra( *fout, specs_crt, crtplots_spill, crtsels_spill, rank, nRanks);

  fout->Write();
  fout->Close();
}