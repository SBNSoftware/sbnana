// Recommend running with a large --stride argument to just use the first file
// in the dataset

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/HistAxis.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"

using namespace ana;

// Since I'm in the SBND group the SAM dataset works fine for me
//const std::string sbnd_wildcard = "official_MCP2021A_CAF_prodgenie_nu_singleinteraction_tpc_sbnd_concat_caf_sbnd";// "official_MCP2021A_CAF_prodoverlay_corsika_cosmics_proton_genie_nu_spill_tpc_sbnd_concat_caf_sbnd";//workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_flat_caf_sbnd";

//const std::string sbnd_wildcard = "/pnfs/icarus/scratch/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_42_00/MultiSigmaAdded/*.root"; // actually icarus NuMI

const std::string sbnd_wildcard = "/pnfs/icarus/scratch/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_37_01_03p02/MultiSigmaAdded_FromConcat/*.root"; // concat icarus files

const std::string state_fname = "test_tracks_state.root";

void test_tracks(bool reload = false)
{
  if(reload || TFile(state_fname.c_str()).IsZombie()){
    SpectrumLoader loader(sbnd_wildcard);

    const TrackVar kLength = SIMPLETRACKVAR(len);
    const TrackCut kLengthCut = kLength>100;

    const Binning binsLength = Binning::Simple(30, 0, 300);

    const TrackHistAxis axLength("RecoLength", binsLength, kLength);

    auto& src = loader.Slices().Tracks();

    Spectrum sAll(src, axLength);
    Spectrum sLong(src[kLengthCut], axLength);

    loader.Go();

    TFile fout(state_fname.c_str(), "RECREATE");
    sAll.SaveTo(&fout, "all");
    sLong.SaveTo(&fout, "long");
  }

  TFile fin(state_fname.c_str());
  Spectrum* sAll = LoadFrom<Spectrum>(&fin, "all").release();
  Spectrum* sLong = LoadFrom<Spectrum>(&fin, "long").release();

  sAll->ToTH1(kPOTnominal, kRed)->Draw("hist");
  sLong->ToTH1(kPOTnominal, kBlue)->Draw("hist same");

  // Redraw the nominals over the top
  sAll->ToTH1(kPOTnominal, kRed)->Draw("hist same");
  sLong->ToTH1(kPOTnominal, kBlue)->Draw("hist same");

  gPad->Print("test_tracks.pdf");
}
