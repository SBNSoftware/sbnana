#include "TFile.h"
#include "TH1.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
// NuMINumuXSec
#include "sbnana/SBNAna/CRTPMTMatching/ICARUSCRTPMTMatching.h"

using namespace ICARUSCRTPMTMatching;

void test_CRTPMT(){

  cpmt.SetGateType(NUMI);

  SpectrumLoader loader("/pnfs/icarus/persistent/users/jskim/data/run_8515/flatcaf/v09_63_00_02/221212_UseOldFlashMatching_FixCAFT1/NUMIMAJORITY/flatcaf_0.root");
  Spectrum *s = new Spectrum("CRTPMTTime", Binning::Simple(3000, -0.15, 0.15), loader, spillvarCRTPMTTime, spillcutCRTPMTCosmicByID);
  loader.Go();

  TFile *f_out = new TFile("output.root","RECREATE");
  f_out->cd();
  TH1D* h = s->ToTH1(6.0e20);
  h->Write();
  f_out->Close();

}
