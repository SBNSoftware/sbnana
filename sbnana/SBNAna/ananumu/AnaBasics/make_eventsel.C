// Make a few spectra with different cuts.

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "cafanacore/Spectrum.h"

#include "eventsel.h"

#include "TFile.h"

#include <iostream>

using namespace ana;

void make_eventsel()
{
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  const std::string fname = "/sbnd/app/users/psihas/SBNCAF-dev/reco-a1e5cc95-e889-474e-87a6-ffc68ec3825f.caf.root";
  //  const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/workshopdemo/numu.caf.root";

  // Source of events
  SpectrumLoader loader(fname);

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum *specs[kNSel][kNVar];

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
      specs[iSel][jVar] = new Spectrum(loader.Slices()[sels[iSel].cut], HistAxis(plots[jVar].label, plots[jVar].bins, plots[jVar].var));
    }
  }

  // This is the call that actually fills in the spectra
  loader.Go();

  // Save spectra to a file, so we can plot them later
  TFile fout("out_demo.root", "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
      std::string mysuffix = sels[iSel].suffix + "_" + plots[jVar].suffix;
      std::cout << "Saving spctra: " << mysuffix << std::endl;
      specs[iSel][jVar]->SaveTo(&fout, mysuffix);
    }
  }

}
