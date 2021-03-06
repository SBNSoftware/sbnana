// Make a few spectra with different cuts.

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "eventsel.h"

#include "TFile.h"

using namespace ana;

void make_eventsel()
{
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  //const std::string fname = "/sbnd/app/users/psihas/SBNCAF-dev/reco-a1e5cc95-e889-474e-87a6-ffc68ec3825f.caf.root";
  //  const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/workshopdemo/numu.caf.root";
  //const std::string fname = "/pnfs/sbnd/persistent/sbndpro/mcp/mc/workshop/SBNWorkshop0421/prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc/v09_19_00_01/caf/flat_caf_3-dac91e87-fdad-4636-855c-107f4a5037d1.root";
  const std::string fname = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_caf_sbnd";

  // Source of events
  SpectrumLoader loader(fname);

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum *specs[kNSel][kNVar];

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
      specs[iSel][jVar] = new Spectrum(plots[jVar].label, plots[jVar].bins, loader, plots[jVar].var, sels[iSel].cut);
    }
  }

  // This is the call that actually fills in the spectra
  loader.Go();

  // Save spectra to a file, so we can plot them later
  //TFile fout("out_demo.root", "RECREATE");
  TFile fout("NOInTime_AllCuts.root", "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
      std::string mysuffix = sels[iSel].suffix + "_" + plots[jVar].suffix;
      std::cout << "Saving spctra: " << mysuffix << std::endl;
      specs[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
    }
  }

}
