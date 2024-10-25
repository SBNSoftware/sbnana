#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"

#include "TFile.h"

using namespace ana;

void MakeTree(){

  // Define inputfiles
  std::vector<std::string> vec_inputs;
  vec_inputs.push_back("/pnfs/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst11Oct24/flatcaf_SBNNuSyst_0.root");

  // SpectrumLoader
  SpectrumLoader loader(vec_inputs);

  // Define Tree variables
  std::vector<std::string> vec_labels_SelectedEvents = {
    "CutType/i",
    "RecoMuonP", "TrueMuonP",
  };
  std::vector<Var> vec_vars_SelectedEvents = {
    kNuMICutType,
    kNuMIMuonCandidateRecoP, kNuMIMuonTrueP,
  };

  Tree *tree_SelectedEvents = new ana::Tree(
    "SelectedEvents",
    vec_labels_SelectedEvents,
    loader,
    vec_vars_SelectedEvents,
    kNoSpillCut,
    kNoCut,
    kNoShift,
    true, true
  );

  loader.Go();

  // Output
  TFile *f_out = new TFile("output.root", "RECREATE");
  tree_SelectedEvents->SaveTo(f_out);
  //tree_SelectedEvents->SaveTo(f_out->GetDirectory(""));

}
