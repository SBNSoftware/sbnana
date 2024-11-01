#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"

#include "TFile.h"

#include "GUNDAMUtils.h"

using namespace ana;

void MakeTree(){

  // Define inputfiles
  std::vector<std::string> vec_inputs;
  vec_inputs.push_back("/pnfs/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst11Oct24/flatcaf_SBNNuSyst_0.root");

  // SpectrumLoader
  SpectrumLoader loader(vec_inputs);

  // Define Tree variables
  // Example here are defined at sbnana/SBNAna/Vars/NuMIXSecVars.h
  std::vector<std::string> vec_labels_SelectedEvents = {
    "CutType/i",
    "RecoMuonP", "TrueMuonP",
    "RecoMuonCos", "TrueMuonCos",
  };
  std::vector<Var> vec_vars_SelectedEvents = {
    kNuMICutType,
    kNuMIMuonCandidateRecoP, kNuMIMuonTrueP,
    kNuMIRecoCosThVtx, kNuMITrueCosThVtx,
  };

  // Event Tree
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
  // NSigmasTree; splines
  // 1) First we need to define ISysts
  std::vector<std::string> NSigmasPsetNames;
  std::vector<const ISyst*> NSigmasISysts;
  std::vector<std::vector<double>> NSigmas;

  std::vector<std::string> GENIEMultisigmaKnobNames = ICARUSGUNDAM::GetGENIEMultisigmaKnobNames();
  std::string SystProviderPrefix = "GENIEReWeight_ICARUS_v2";
  for(unsigned int i=0; i<GENIEMultisigmaKnobNames.size(); i++){
    NSigmasPsetNames.push_back( GENIEMultisigmaKnobNames.at(i) );

    std::string psetname = SystProviderPrefix+"_multisigma_"+GENIEMultisigmaKnobNames.at(i);
    NSigmasISysts.push_back( new SBNWeightSyst(psetname) );

    NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
  }

  NSigmasTree *tree_SelectedEvents_NSigmas = new ana::NSigmasTree(
    "SelectedEvents_NSigmas",
    NSigmasPsetNames,
    loader,
    NSigmasISysts,
    NSigmas,
    kNoSpillCut,
    kNoCut,
    kNoShift,
    true, true
  );

  loader.Go();

  // Output
  TFile *f_out = new TFile("output.root", "RECREATE");
  tree_SelectedEvents->SaveTo(f_out);
  tree_SelectedEvents_NSigmas->SaveTo(f_out);

}
