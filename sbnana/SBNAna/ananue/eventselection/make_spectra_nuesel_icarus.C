//! ////////////////////////////////////////////////////////////////////////////
//! Authors: Diana Patricia Mendez, Jacob Smith (smithja)                     //
//! Contact: jacob.a.smith@stonybrook.edu                                     //
//! Last edited: February 24th, 2025                                          //   
//!                                                                           // 
//! Makes Spectra of the variables defined in helper with specific            // 
//! cuts applied to the sample(s).                                            //
//! Takes a CAF and outputs a ROOT file containing the resulting histograms.  //
//! ////////////////////////////////////////////////////////////////////////////

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "helper_nuesel_icarus.h"
#include "tools_nuesel_icarus.h"

#include "TFile.h"
#include "TTreeReader.h"

using namespace ana;

void make_spectra_nuesel_icarus( std::string finname = "nus", 
                                 int setno = 1, 
                                 bool selspill = false) {
  std::string settag = std::to_string(setno);

  //! Create input and output directories and files.
  std::string findir = "/inputdir_path/";
  std::string *pindir = &findir;
  create_dir( pindir);
  std::string foutdir = "outputs/";
  std::string *poutdir = &foutdir;
  create_dir( poutdir);

  const std::string finsample = "refactored_g4_icarus_step2_g4step2_detsim_"
                                "stage0_stage1_13629324_0.flat.caf"
                                "-e8a8c7bb-7ad3-46b2-a5a6-2496ffc7199f.root"; 
//  const std::string finsample = findir + finname + "_hadded" + settag + ".flatcaf.root";
//  const std::string finsample = "icaruspro_production_v09_89_01_01_2024A_ICARUS_MC_Sys_NuCos_2024A_MC_Sys_NuCos_CV_flatcaf";

//! Create the output filename with the "spill" tag if setno is declared.
//! @note Since setno has a default value of 1 in the make_spectra_nuesel_icarus()
//!       function, you will likely have to delete setno to use set the output
//!       file name to have the "slice" tag.
  const std::string foutname = foutdir + finsample + "_OUTPUT_SPECTRA.root";
//  const std::string foutname = foutdir + finname + "_spectra_hadded" + settag + (selspill ? "_spill" : "_slice") + ".root";

  SpectrumLoader loader( finsample);

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

  std::vector<PlotDef> plots;
  std::vector<SelDef> types;
  std::vector<SelDef> sels;

  //! Select all slices or a single slice per spill.
  if( selspill){
    //! @note may not work due to differences in structure of slice/spill objects
    plots = plots_spill;
    types = types_spill;
    sels  = sels_spill;
  }
  else{ // using slices
    plots = plots_slice;
    types = types_slice;
    sels  = sels_slice;
  }

  //! Since we're making Spectra of all the different combinations of variables, 
  //! selection criteria, and interaction type, define the total numbers of
  //! these Spectra:
  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  const unsigned int kNType = types.size();

  //! @note These variables are declared in helper_nuesel_icarus.h
  const unsigned int kNVarCRTSpill = crtplots_spill.size();
  const unsigned int kNSelCRTSpill = crtsels_spill.size();

  //! Make an array of Spectra with each selection-type-var combination:

  //! @note Declaring these Spectra as arrays of pointers likely improves
  //!       the speed of this program by A LOT given the number of Spectra
  //!       we're making.
  //!
  //!       HOWEVER, doing so leads to memory leaks. It does not seem to affect
  //!       the output file containing all of the Spectra. The error, a 
  //!       "double free or corruption (!prev)" error, indicates this.
  //!       
  //!       I (Jacob Smith, smithja) RECOMMEND TROUBLESHOOTING THESE MEMORY
  //!       ERRORS BEFORE INCORPORATING THIS SCRIPT INTO A LARGER FRAMEWORK.
  //!       This would likely be done easiest with std::vectors.

  Spectrum *specs[kNSel][kNType][kNVar];
  Spectrum *specsveto[kNSel][kNType][kNVar];

  Spectrum *specscrtspill[kNSelCRTSpill][kNVarCRTSpill];


  for( unsigned int iSel = 0; iSel < kNSel; ++iSel) {
    for( unsigned int iType = 0; iType < kNType; ++iType) {
      for( unsigned int iVar = 0; iVar < kNVar; ++iVar) {
        //! Spectra WITHOUT CRT veto:
        specs[iSel][iType][iVar] = new Spectrum(
          plots[iVar].label, 
          plots[iVar].bins, 
          loader, 
          plots[iVar].var, 
          sels[iSel].cut && types[iType].cut );

        //! Spectra WITH CRT veto (kCRTHitVetoFD):
        specsveto[iSel][iType][iVar] = new Spectrum(
          plots[iVar].label, 
          plots[iVar].bins, 
          loader, 
          plots[iVar].var, 
          kCRTHitVetoFD, 
          sels[iSel].cut && types[iType].cut );
      }
    }
  }

  for( unsigned int iSel = 0; iSel < kNSelCRTSpill; ++iSel) {
    for( unsigned int iVar = 0; iVar < kNVarCRTSpill; ++iVar) {
      specscrtspill[iSel][iVar] = new Spectrum(
        crtplots_spill[iVar].label, 
        crtplots_spill[iVar].bins, 
        loader, 
        crtplots_spill[iVar].var, 
        crtsels_spill[iSel].cut);
    }
  }

  //! This line is what actually fills all of the Spectra.
  loader.Go();

  //! Save the Spectra to a file so we can plot them later.
  TFile fout(foutname.c_str(), "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for(unsigned int iType = 0; iType < kNType; ++iType){
      for( unsigned int iVar = 0; iVar < kNVar; ++iVar ){
        Spectrum *spec      = specs[iSel][iType][iVar];
        Spectrum *spec_veto = specsveto[iSel][iType][iVar];

        std::string mysuffix = types[iType].suffix + "_" + 
                               sels[iSel].suffix + "_" + 
                               plots[iVar].suffix;
        std::string mysuffixveto = types[iType].suffix + "_" + 
                                   sels[iSel].suffix + "_veto_" + 
                                   plots[iVar].suffix;

        std::cout << "Saving spectra: " << mysuffix << std::endl;

        ////////////////////////////////////////////////////////////////////////
        //! Pseudo-scale Cosmics
        
        //! If we are dealing with cosmic data, apply cosmic POT defined above.
        if( finname == "cosmics") {
          std::cout << "Fake POT scale cosmics!" << std::endl;
          const double thisPOT( spec->POT());

          //! If this Spectra's POT is smaller than the minimum perceptible
          //! difference for C++, overwrite the POT as defined above.
          if( thisPOT < std::numeric_limits<double>::epsilon()){
            spec->OverridePOT( cosmicPOT);
            spec_veto->OverridePOT( cosmicPOT);
          }
	      }
        ////////////////////////////////////////////////////////////////////////

        spec->SaveTo( fout.mkdir( mysuffix.c_str()));
        spec_veto->SaveTo( fout.mkdir( mysuffixveto.c_str()));
      }
    }
  }

  for( unsigned int iSel = 0; iSel < kNSelCRTSpill; ++iSel ){
    for( unsigned int iVar = 0; iVar < kNVarCRTSpill; ++iVar ){
      std::string mysuffix = crtsels_spill[iSel].suffix + "_" + 
                             crtplots_spill[iVar].suffix;

      std::cout << "Saving spectra: " << mysuffix << std::endl;

      specscrtspill[iSel][iVar]->SaveTo( fout.mkdir( mysuffix.c_str()));

    }
  }
  fout.Close();
}
