//////////////////////////////////////////////////////////////////////////////
// Authors: Diana Patricia Mendez, Jacob Smith (smithja)                    //
// Contact: jacob.a.smith@stonybrook.edu                                    //
// Last edited: February 24th, 2025                                         //   
//                                                                          // 
// Makes Spectra of the variables defined in helper with specific           // 
// cuts applied to the sample(s).                                           //
// Takes a CAF and outputs a ROOT file containing the resulting histograms. //
//////////////////////////////////////////////////////////////////////////////

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

//  std::string findir = "/inputdir_path/";
//  std::string *pindir = &findir;
//  create_dir( pindir);

  std::string foutdir = "outputs/";
//  std::string *poutdir = &foutdir;
//  create_dir( poutdir);

//  const std::string finsample = findir + finname + 
//                                "_hadded" + settag + ".flatcaf.root";
//  const std::string finsample = "icaruspro_production_v09_89_01_01_2024A_"
//                                "ICARUS_MC_Sys_NuCos_2024A_MC_Sys_NuCos_CV_flatcaf";
//  const std::string foutname = foutdir + finname + 
//                               "_spectra_hadded" + settag + 
//                               (selspill ? "_spill" : "_slice") + ".root";

  const std::string finsample = "refactored_g4_icarus_step2_g4step2_detsim_stage0_stage1_13629324_0.flat.caf-e8a8c7bb-7ad3-46b2-a5a6-2496ffc7199f.root"; 
  const std::string foutname = foutdir + "TEST_OUTPUT.root";

  // Source of events from a flatcaf:
  SpectrumLoader loader( finsample);

/*
  ////////////////////////////////
  //// pseudo scale cosmics

  TFile f( finsample.c_str());
  TTree* myTree = (TTree*) f.Get("recTree");
  TTreeReader cafReader("recTree", &f);

//  TTreeReaderValue< unsigned int> run(     cafReader, "hdr.run");
//  TTreeReaderValue< unsigned int> subrun(  cafReader, "hdr.subrun");
//  TTreeReaderValue< unsigned int> ngenevt( cafReader, "hdr.ngenevt");

  // Use this instead of the above if using flatcafs.
  // Otherwise cosmics will always have 'POT' set to zero.
  TTreeReaderValue< unsigned int> run(     cafReader, "rec.hdr.run");
  TTreeReaderValue< unsigned int> subrun(  cafReader, "rec.hdr.subrun");
  TTreeReaderValue< unsigned int> ngenevt( cafReader, "rec.hdr.ngenevt");

  unsigned int totalGenEvt(0);

  std::map<std::pair< unsigned int, unsigned int>, unsigned int> runSubRunMap;
  while( cafReader.Next()) {
    std::pair< unsigned int, unsigned int> thisPair(*run, *subrun);
    if( runSubRunMap[thisPair] == 0) {
      ++runSubRunMap[thisPair];
      totalGenEvt += *ngenevt;
    }
  }

  const unsigned int numTriggers( myTree->GetEntries());

  std::cout << "Num Triggers: " << numTriggers << std::endl;
  std::cout << "Num Triggers from generated events: " << totalGenEvt << std::endl;

  const float nuPerSpill( 1.f / 21.8);
  const float POTperSpill( 5e12);
  const float cosmicPOT( totalGenEvt * POTperSpill / (1.f - nuPerSpill));

  std::cout << "cosmicPOT: " << cosmicPOT << std::endl;
  const float spills_per_cos_trig = (cosmicPOT / numTriggers) / POTperSpill; 
  std::cout << "Spills per cosmic Trigger: " << spills_per_cos_trig << std::endl;

  //// end pseudo scale cosmics
  ////////////////////////////////
*/

  // select all slices or a single slice per spill
  std::vector<PlotDef> plots;
  std::vector<SelDef> types;
  std::vector<SelDef> sels;

//  if( selspill){
    // may not work due to differences in structures of slice and spill objects
//    plots = plots_spill;
//    types = types_spill;
//    sels  = sels_spill
//  }
//  else{ // using slices
    plots = plots_slice;
    types = types_slice;
    sels  = sels_slice;
//  }

  // Since we're making Spectra of all the different combinations of variable, 
  // selection criteria, and interaction type, define the total numbers of
  // these Spectra:
  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  const unsigned int kNType = types.size();
  const unsigned int kNVarCRTSpill = crtplots_spill.size();
  const unsigned int kNSelCRTSpill = crtsels_spill.size();

  // Make an array of Spectra with each selection-type-var combination:

  // NOTE: Declaring these Spectra as arrays of pointers likely improves
  //       the speed of this program by A LOT given the number of Spectra
  //       we're making.
  //
  //       HOWEVER, doing so leads to memory leaks. It does not seem to affect
  //       the output file containing all of the Spectra. The error, a 
  //       "double free or corruption (!prev)" error, indicates this.
  //       
  //       I (Jacob Smith, smithja) RECOMMEND TROUBLESHOOTING THESE MEMORY
  //       ERRORS BEFORE INCORPORATING THIS SCRIPT INTO A LARGER FRAMEWORK.
  //       This would likely be done easiest with std::vectors.

  Spectrum *specs[kNSel][kNType][kNVar];
  Spectrum *specsveto[kNSel][kNType][kNVar];
  Spectrum *specscrtspill[kNSelCRTSpill][kNVarCRTSpill];

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel) {
    for( unsigned int iType = 0; iType < kNType; ++iType) {
      for( unsigned int jVar = 0; jVar < kNVar; ++jVar) {
        // without CRT veto:
        specs[iSel][iType][jVar] = new Spectrum(
          plots[jVar].label, 
          plots[jVar].bins, 
          loader, 
          plots[jVar].var, 
          sels[iSel].cut && types[iType].cut );

        // with CRT veto:
        specsveto[iSel][iType][jVar] = new Spectrum(
          plots[jVar].label, 
          plots[jVar].bins, 
          loader, 
          plots[jVar].var, 
          kCRTHitVetoFD, 
          sels[iSel].cut && types[iType].cut );
      }
    }
  }

  for( unsigned int iSel = 0; iSel < kNSelCRTSpill; ++iSel) {
    for( unsigned int jVar = 0; jVar < kNVarCRTSpill; ++jVar) {
      specscrtspill[iSel][jVar] = new Spectrum(
        crtplots_spill[jVar].label, 
        crtplots_spill[jVar].bins, 
        loader, 
        crtplots_spill[jVar].var, 
        crtsels_spill[iSel].cut);
    }
  }
  // This line is what actually fills all of the Spectra:
  loader.Go();

  // Save the Spectra to a file so we can plot them later:
  TFile fout(foutname.c_str(), "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for(unsigned int iType = 0; iType < kNType; ++iType){
      for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
        std::string mysuffix = types[iType].suffix + "_" + 
                               sels[iSel].suffix + "_" + 
                               plots[jVar].suffix;
        std::string mysuffixveto = types[iType].suffix + "_" + 
                                   sels[iSel].suffix + "_veto_" + 
                                   plots[jVar].suffix;

        std::cout << "Saving spectra: " << mysuffix << std::endl;

/*

        ////////////////////////////////
        //// pseudo scale cosmics

        // If we are dealing with cosmic data, apply cosmic POT defined above.
        if( finname == "cosmics") {
          std::cout << "Fake POT scale cosmics!" << std::endl;
          const double thisPOT( spec->POT());

          // If this Spectra's POT is smaller than the minimum perceptible
          // difference for C++, overwrite the POT as defined above.
          if( thisPOT < std::numeric_limits<double>::epsilon()){
            spec->OverridePOT( cosmicPOT);
            specveto->OverridePOT( cosmicPOT);
          }
	    } // fake cosmics pot

        //// end pseudo scale cosmics
        ////////////////////////////////

*/ 
        specs[iSel][iType][jVar]->SaveTo( fout.mkdir( mysuffix.c_str()));
        specsveto[iSel][iType][jVar]->SaveTo( fout.mkdir( mysuffixveto.c_str()));

//        delete specs[iSel][iType][jVar];
//        delete specsveto[iSel][iType][jVar];
      }
    }
  }


  for( unsigned int iSel = 0; iSel < kNSelCRTSpill; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVarCRTSpill; ++jVar ){
      std::string mysuffix = crtsels_spill[iSel].suffix + "_" + 
                             crtplots_spill[jVar].suffix;

      std::cout << "Saving spectra: " << mysuffix << std::endl;

      specscrtspill[iSel][jVar]->SaveTo( fout.mkdir( mysuffix.c_str()));

//      delete specscrtspill[iSel][jVar];
    }
  }

  fout.Close();
//  delete &fout;
}
