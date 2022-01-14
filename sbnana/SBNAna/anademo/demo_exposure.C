#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "sbnana/SBNAna/Vars/NumuVarsSBND202106.h"

#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"

using namespace ana;

// Location of intime cosmics files
const std::string cosmic_wildcard = "workshop_SBNWorkshop0421_prodcorsika_proton_intime_filter_flat_caf_sbnd";

// Location of beam (plus out of time cosmics overlay) files
const std::string beam_wildcard = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_flat_caf_sbnd";

// For SAM to find these files, you need to be in the SBND group

void demo_exposure()
{
  SpectrumLoader loaderBeam(beam_wildcard);
  SpectrumLoader loaderCosmic(cosmic_wildcard);

  const Cut kNumuSel = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut;

  const Cut kIsCosmic = SIMPLEVAR(truth.index) < 0;

  const Binning binsLen = Binning::Simple(30, 0, 600);

  const HistAxis axLen("Primary muon track length (cm)", binsLen, kPrimaryMuonTrkLen);

  // Beam components
  Spectrum sBeamCC(loaderBeam, axLen, kNumuSel && kIsNumuCC);
  Spectrum sBeamNC(loaderBeam, axLen, kNumuSel && kIsNC);
  // Cosmics in the beam readouts ("out of time cosmics")
  Spectrum sBeamCosmic(loaderBeam, axLen, kNumuSel && kIsCosmic);
  // Cosmics outside beam readouts ("in time cosmics")
  Spectrum sCosmic(loaderCosmic, axLen, kNumuSel);

  loaderBeam.Go();
  loaderCosmic.Go();

  // These three can be plotted correctly just by scaling to a target POT
  TH1* hBeamCC = sBeamCC.ToTH1(kPOTnominal);
  TH1* hBeamNC = sBeamNC.ToTH1(kPOTnominal, kBlue);
  TH1* hBeamCosmic = sBeamCosmic.ToTH1(kPOTnominal, kOrange);

  // For the in-time cosmics:
  //
  // Figure out how much livetime (number of spills) this POT corresponds to
  const double nominalIntensity = 5e12; // POT per spill
  const double nominalLivetime = kPOTnominal/nominalIntensity;
  // Find out how many spills are already accounted for by the beam MC. This is
  // equivalent to sBeamCC.Livetime()*kPOTnominal/sBeamCC.POT()
  const double beamLivetime = sBeamCC.FakeData(kPOTnominal).Livetime();
  // How many spills still need to be accounted for out of the cosmics files
  const double cosmicLivetime = nominalLivetime - beamLivetime;

  TH1* hCosmic = sCosmic.ToTH1(cosmicLivetime, kMagenta, kSolid,
                               kLivetime /* treat the first argument as a livetime */);

  hBeamCC->Draw("hist");
  hBeamNC->Draw("hist same");
  hBeamCosmic->Draw("hist same");
  hCosmic->Draw("hist same");

  TLegend* leg = new TLegend(.5, .5, .85, .85);
  leg->SetFillStyle(0);
  leg->AddEntry(hBeamCC, "#nu_{#mu} CC", "l");
  leg->AddEntry(hBeamNC, "NC", "l");
  leg->AddEntry(hBeamCosmic, "Out-of-time cosmics", "l");
  leg->AddEntry(hCosmic, "In-time cosmics", "l");
  leg->Draw();

  std::cout << "Scaling to " << kPOTnominal << " POT and "
            << nominalLivetime << " spills (assuming "
            << nominalIntensity << " POT/spill):" << std::endl;
  std::cout << "  Numu CCs           : " << sBeamCC.Integral(kPOTnominal) << std::endl;
  std::cout << "       NCs           : " << sBeamNC.Integral(kPOTnominal) << std::endl;
  std::cout << "  out-of-time cosmics: " << sBeamCosmic.Integral(kPOTnominal) << std::endl;
  std::cout << "      in-time cosmics: " << sCosmic.Integral(cosmicLivetime, 0, kLivetime) << std::endl;

  gPad->Print("demo_exposure.png");
}
