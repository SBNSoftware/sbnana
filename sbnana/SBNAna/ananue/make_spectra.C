// Make a few spectra with different cuts.

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "cafanacore/Spectrum.h"

#include "helper.h"

#include "TFile.h"

#include <iostream>

using namespace ana;

void make_spectra()
{
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/Nov2020CAFs/hadded_1_nuecosmics.caf.root";

  // Source of events
  SpectrumLoader loader(fname);

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  const unsigned int kNVarSpill = plots_spill.size();
  const unsigned int kNSelSpill = sels_spill.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum *specs[kNSel][kNVar]; // (Slice)Var with (Slice)Cut

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
      specs[iSel][jVar] = new Spectrum(loader.Slices()[sels[iSel].cut], HistAxis(plots[jVar].label, plots[jVar].bins, plots[jVar].var));
    }
  }

  // Example spectra with both SpillCut and Cut
  //  Spectrum *s0 = new Spectrum(loader[kFlashTrigger], SpillHistAxis("crthitx__flashtrig", kPositionXFDBinning, kCRTHitX)); // SpillMultiVar with SpillCut
  Spectrum *s1 = new Spectrum(loader[kFirstEvents].Slices()[!kShortShower], HistAxis("shwlen__firstevt_longshw", kLengthBinning, kRecoShower_Length)); // (Slice)Var with SpillCut and (Slice)Cut
  Spectrum *s2 = new Spectrum(loader[kFlashTrigger].Slices()[!kShortShower], HistAxis("shwlen__flashtrig_longshw", kLengthBinning, kRecoShower_Length)); // (Slice)Var with SpillCut and (Slice)Cut
  Spectrum *s3 = new Spectrum(loader[kFlashTrigger].Slices()[kShortShower], HistAxis("shwlen__flashtrig_shortshw", kLengthBinning, kRecoShower_Length)); // (Slice)Var with SpillCut and (Slice)Cut

  // This is the call that actually fills in the spectra
  loader.Go();

  // Save spectra to a file, so we can plot them later
  TFile fout("spectra_eventsel_nue.root", "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
      std::string mysuffix = sels[iSel].suffix + "_" + plots[jVar].suffix;
      std::cout << "Saving spectra: " << mysuffix << std::endl;
      specs[iSel][jVar]->SaveTo(&fout, mysuffix);
    }
  }

  //  s0->SaveTo(&fout, "s0");
  s1->SaveTo(&fout, "s1");
  s2->SaveTo(&fout, "s2");
  s3->SaveTo(&fout, "s3");

}
