#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "TH1.h"
#include "TPad.h"

using namespace ana;

// Since I'm in the SBND group the SAM dataset works fine for me
const std::string sbnd_wildcard = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_flat_caf_sbnd";

void spectra_sbnd()
{
  SpectrumLoader loader(sbnd_wildcard);

  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kRecoE = kTrueE; // TODO

  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSel = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut;

  const Cut kIsCosmic = SIMPLEVAR(truth.index) < 0;

  // This looks much better than the official fit binning 
  const Binning binsEnergy = Binning::Simple(30, 0, 3);

  // TODO should make this an official function somewhere
  /*
  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                                   0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                   1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);
  */

  const HistAxis axEnergy(/*Reconstructed*/"True energy (GeV)", binsEnergy, kRecoE);

  Spectrum sTot(loader, axEnergy, kNumuSpillSel, kNumuSel);
  Spectrum sNC(loader, axEnergy, kNumuSpillSel, kNumuSel && kIsNC);
  //  Spectrum sCosmic(loader, axEnergy, kNumuSpillSel, kNumuSel && kIsCosmic);

  loader.Go();

  sTot.ToTH1(kPOTnominal)->Draw("hist");
  TH1* hNC = sNC.ToTH1(kPOTnominal, kBlue);
  hNC->SetFillColor(kBlue-10);
  hNC->Draw("histf same");
  //  TH1* hCosmic = sCosmic.ToTH1(kPOTnominal, kMagenta);
  //  hCosmic->Draw("hist same");

  std::cout << sTot.Integral(kPOTnominal) << " events total" << std::endl;
  std::cout << "  of which " << sNC.Integral(kPOTnominal) << " NC" << std::endl;
  //  std::cout << "  and " << sCosmic.Integral(kPOTnominal) << " cosmics" << std::endl;

  gPad->Print("spectrum_sbnd.pdf");
}
