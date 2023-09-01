// Make a few spectra with different cuts.

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "selection_helper.h"

#include "TFile.h"

using namespace ana;

void make_spectra()
{
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  //const std::string fname = "/icarus/data/users/rhowell/nucosmics_hadded.root";
  const std::string fname = "./pfpnew_tree.flat.root";

  // Source of events
  SpectrumLoader loader(fname);

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  
  const unsigned int kNSpillVar = spillPlots.size();
  const unsigned int kNSpillSel = spillSels.size();

  const unsigned int kNMultiVars = multiPlots.size();
  const unsigned int kNMultiSels = sels.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum *specs[kNSel][kNVar];

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
      specs[iSel][jVar] = new Spectrum(plots[jVar].label, plots[jVar].bins, loader, plots[jVar].var, sels[iSel].cut);
    }
  }

  Spectrum *spillSpecs[kNSpillSel][kNSpillVar];

  for(unsigned int iSel = 0; iSel < kNSpillSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNSpillVar; ++jVar){
      spillSpecs[iSel][jVar] = new Spectrum(spillPlots[jVar].label, spillPlots[jVar].bins, loader, spillPlots[jVar].var, spillSels[iSel].cut, kSpillUnweighted);
    }
  }

  Spectrum *multiSpecs[kNMultiSels][kNMultiVars];

  for(unsigned int iSel = 0; iSel < kNMultiSels; ++iSel){
    for(unsigned int jVar = 0; jVar < kNMultiVars; ++jVar){
      multiSpecs[iSel][jVar] = new Spectrum(multiPlots[jVar].Xlabel, multiPlots[jVar].Ylabel, loader, multiPlots[jVar].binsX, multiPlots[jVar].varX, multiPlots[jVar].binsY, multiPlots[jVar].varY, kNoSpillCut,sels[iSel].cut);
    }
  }

  // This is the call that actually fills in the spectra
  loader.Go();

  // Save spectra to a file, so we can plot them later
  TFile fout("nue_cut_distros.root", "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
      std::string mysuffix = sels[iSel].suffix + "_" + plots[jVar].suffix;
      std::cout << "Saving spectra: " << mysuffix << std::endl;
      specs[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
    }
  }

  for( unsigned int iSel = 0; iSel < kNSpillSel; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNSpillVar; ++jVar ){
      std::string mysuffix = spillSels[iSel].suffix + "_" + spillPlots[jVar].suffix;
      std::cout << "Saving spectra: " << mysuffix << std::endl;
      spillSpecs[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
    }
  }

  for( unsigned int iSel = 0; iSel < kNMultiSels; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNMultiVars; ++jVar ){
      std::string mysuffix = sels[iSel].suffix + "_" + multiPlots[jVar].suffix;
      std::cout << "Saving multi spectra: " << mysuffix << std::endl;
      multiSpecs[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
    }
  }
}

