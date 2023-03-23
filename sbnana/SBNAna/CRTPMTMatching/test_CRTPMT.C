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

  //SpectrumLoader loader("/pnfs/icarus/persistent/users/jskim/data/run_8515/flatcaf/v09_63_00_02/221212_UseOldFlashMatching_FixCAFT1/NUMIMAJORITY/flatcaf_0.root");
  SpectrumLoader loader("/pnfs/icarus/persistent/users/jskim/Run2/mc/NUMI_Nu_Cosmics/flatcaf/v09_63_00_02/230114_AddGHEPPtl_FixCRTTimingInCAF_CRTPMTMatchingVar/G18_10a_02_11a/flatcaf_0.root");

  // if MC,
  cpmt.UseTS0 = true;

  Spectrum *s_CRTPMTTime = new Spectrum("CRTPMTTime", Binning::Simple(3000, -0.15, 0.15), loader, spillvarCRTPMTTime, kNoSpillCut);
  Spectrum *s_CRTPMTMatchingID = new Spectrum("CRTPMTMatchingID", Binning::Simple(15, 0., 15.), loader, spillvarCRTPMTMatchingID, kNoSpillCut);
  Spectrum *s_MatchID2_CRTHitPosXs = new Spectrum("MatchID2_CRTHitPosXs", Binning::Simple(800, -400, 400.), loader, spillvarMatchID2_CRTHitPosXs, kNoSpillCut);

  loader.Go();

  TFile *f_out = new TFile("output.root","RECREATE");
  f_out->cd();

  TH1D* h_CRTPMTTime = s_CRTPMTTime->ToTH1(6.0e20);
  h_CRTPMTTime->SetName("CRTPMTTime");
  h_CRTPMTTime->Write();

  TH1D* h_CRTPMTMatchingID = s_CRTPMTMatchingID->ToTH1(6.0e20);
  h_CRTPMTMatchingID->SetName("CRTPMTMatchingID");
  h_CRTPMTMatchingID->Write();

  TH1D* h_MatchID2_CRTHitPosXs = s_MatchID2_CRTHitPosXs->ToTH1(6.0e20);
  h_MatchID2_CRTHitPosXs->SetName("MatchID2_CRTHitPosXs");
  h_MatchID2_CRTHitPosXs->Write();

  f_out->Close();

}
