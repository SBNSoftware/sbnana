#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"

#include "TFile.h"
#include "GUNDAMUtils.h"

using namespace ana;

const Cut IsForTree([](const caf::SRSliceProxy* slc) {
  int cuttpye = kNuMICutType(slc);
  if(cuttpye==0) return false;
  else return true;
});

void MakeTree(){

  // Define inputfiles
  std::vector<std::string> vec_inputs;
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/000d2fcf-6e3b-4f86-a84b-ee676dbb89a7-flatcaf_SBNNuSyst_383.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/002799ea-7b9a-4687-a804-91f0bda84090-flatcaf_SBNNuSyst_41.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/003b40f8-1a5c-4f6c-a543-66077c0c8446-flatcaf_SBNNuSyst_29.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/00bc1b3d-5d5b-41cf-b39f-660c355dcee2-flatcaf_SBNNuSyst_427.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/01eb6788-e9e3-4316-a423-2dac41cab8ef-flatcaf_SBNNuSyst_76.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/039f09a1-1551-4042-aaea-eae8be942899-flatcaf_SBNNuSyst_110.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/03bd5678-0472-4e77-8586-e5804fe8b0cd-flatcaf_SBNNuSyst_218.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/04155113-0115-4d87-ba79-a30cd26072ab-flatcaf_SBNNuSyst_345.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/048f63e1-cd90-4981-acfc-3a75b26d028d-flatcaf_SBNNuSyst_86.root");
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/2023A_NuMI/v09_72_00_08_AddingG4Reweight/MC/NuMI_MC_Nu_Phase2_v09_72_00_03p01_deduped_Merged_flatcaf_ReNuSyst250403_WithTrackSplit_BNBFixedProb/053cdbd0-51cf-4cac-8db5-ce86907799e3-flatcaf_SBNNuSyst_464.root");

  bool IsData = false;

  std::string samdef_MC = "jskim_2023A_NuMI_Reproc_v09_72_00_08_CustomG4RW_NuMI_MC_Nu_Phase2_2023ANuMIReproc_CAFTypeCommonRemerge_ReNuSyst250403_WithTrackSplit_BNBFixedProb";

  // SpectrumLoader
  SpectrumLoader loader(vec_inputs); //, kBeam, 10);
  //SpectrumLoader loader(samdef_MC);

  //===================================
  // EVENT SELECTION TREE

  // Define Tree variables
  // Example here are defined at sbnana/SBNAna/Vars/NuMIXSecVars.h
  std::vector<std::string> vec_labels_SelectedEvents = {
    "CutType/i",
    "IsSignal/i",
    "RecoMuonP", "TrueMuonP",
    "RecoMuonCos", "TrueMuonCos",
    "LeadingChargedPionCandidateLength",
    "NuMIFluxCorrection",
  };
  std::vector<Var> vec_vars_SelectedEvents = {
    kNuMICutType,
    kNuMISliceSignalType,
    kNuMIMuonCandidateRecoP, kNuMIMuonTrueP,
    kNuMIRecoCosThVtx, kNuMITrueCosThVtx,
    kNuMILeadingChargedPionCandidateLength,
    kGetNuMIFluxCorrection,
  };

  if(IsData){

    vec_labels_SelectedEvents.push_back( "IsData/i" );
    vec_vars_SelectedEvents.push_back(
      Var([](const caf::SRSliceProxy* slc) -> int {
        return 1;
      })
    );

  }
  else{

    vec_labels_SelectedEvents.push_back( "IsData/i" );
    vec_vars_SelectedEvents.push_back(
      Var([](const caf::SRSliceProxy* slc) -> int {
        return 0;
      })
    );

  }

  // Event Tree

  Tree *tree_SelectedEvents = new ana::Tree(
    "SelectedEvents",
    vec_labels_SelectedEvents,
    loader,
    vec_vars_SelectedEvents,
    kNoSpillCut,
    IsForTree,
    kNoShift,
    true, true
  );

  // Systematics

  std::string SystProviderPrefix = "GENIEReWeight_ICARUS_v2";

  // NSigmasTree; splines

  std::vector<std::string> NSigmasPsetNames;
  std::vector<const ISyst*> NSigmasISysts;
  std::vector<std::vector<double>> NSigmas;

  // - Multisigma (defined for sigma = -3~+3)

  std::vector<std::string> GENIEMultisigmaKnobNames = ICARUSGUNDAM::GetGENIEMultisigmaKnobNames();
  for(unsigned int i=0; i<GENIEMultisigmaKnobNames.size(); i++){
    NSigmasPsetNames.push_back( GENIEMultisigmaKnobNames.at(i) );

    std::string psetname = SystProviderPrefix+"_multisigma_"+GENIEMultisigmaKnobNames.at(i);
    NSigmasISysts.push_back( new SBNWeightSyst(psetname) );

    NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
  }
  // - Morph (defeind between sigma = 0~1; we mirror them to negative values)
  // - Also it was noted that adding sigma = +-0.5 helps making the spline smoother

  std::vector<std::string> GENIEMorphKnobNames = ICARUSGUNDAM::GetGENIEMorphKnobNames();
  for(unsigned int i=0; i<GENIEMorphKnobNames.size(); i++){
    NSigmasPsetNames.push_back( GENIEMorphKnobNames.at(i) );

    std::string psetname = SystProviderPrefix+"_multisigma_"+GENIEMorphKnobNames.at(i);
    NSigmasISysts.push_back( new SBNWeightMirrorSyst(psetname) );

    NSigmas.push_back( {-1, -0.5, 0, 0.5, 1} );
  }
  // - From NuSyst (you may not have this if your CAFs did not run on nusystemaitcs (sbnnusyst)


  std::vector<std::string> GENIENuSystKnobNames = ICARUSGUNDAM::GetNuSystZExpPCANames();
  for(unsigned int i=0; i<GENIENuSystKnobNames.size(); i++){
    NSigmasPsetNames.push_back( GENIENuSystKnobNames.at(i) );

    std::string psetname = "ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_"+GENIEMorphKnobNames.at(i);
    NSigmasISysts.push_back( new SBNWeightSyst(psetname) );

    NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );

  }

  // - Define NSigmasTree

  NSigmasTree *tree_SelectedEvents_NSigmas = new ana::NSigmasTree(
    "SelectedEvents_NSigmas",
    NSigmasPsetNames,
    loader,
    NSigmasISysts,
    NSigmas,
    kNoSpillCut,
    IsForTree,
    kNoShift,
    true, true
  );

  // NUniversesTree

  std::vector<std::string> NUniversesPsetNames;
  std::vector<std::vector<Var>> NUniversesVarVectors;
  std::vector<unsigned int> NUniversesNUnivs;

  std::vector<std::string> GENIEDependentKnobNames = ICARUSGUNDAM::GetGENIEDependentKnobNames();
  std::map<std::string, std::vector<Var>> map_DepDialName_to_UniverseWeights;
  for(const std::string& name: GENIEDependentKnobNames){
    map_DepDialName_to_UniverseWeights[name] = {};
    std::string psetname = SystProviderPrefix+"_multisim_"+name;
    for(int u=0; u<100; u++){
      map_DepDialName_to_UniverseWeights[name].push_back( GetUniverseWeight(psetname, u) );
    }
  }

  for(const std::string& name: GENIEDependentKnobNames){
    NUniversesPsetNames.push_back( name );
    NUniversesVarVectors.push_back( map_DepDialName_to_UniverseWeights[name] );
    NUniversesNUnivs.push_back( 100 );
  }

  // - Define NUniversesTree

  NUniversesTree *tree_SelectedEvents_NUniverses = new ana::NUniversesTree(
    "SelectedEvents_NUniverses",
    NUniversesPsetNames,
    loader,
    NUniversesVarVectors,
    NUniversesNUnivs,
    kNoSpillCut,
    IsForTree,
    kNoShift,
    true, true
  );

  //===================================
  // TRUTH TREE TREE

  std::vector<std::string> vec_labels_Truth = {
    "TrueMuonP",
    "TrueMuonCos",
  };
  std::vector<TruthVar> vec_vars_Truth = {
    kTruth_MuonP,
    kTruth_MuonNuCosineTheta,
  };

  Tree *tree_Truth = new ana::Tree(
    "trueEvents",
    vec_labels_Truth,
    loader,
    vec_vars_Truth,
    kNoSpillCut,
    kTruthCut_IsSignal,
    kNoCut,
    kNoShift,
    true
  );

  // - Define NSigmasTree

  NSigmasTree *tree_Truth_NSigmas = new ana::NSigmasTree(
    "trueEvents_NSigmas",
    NSigmasPsetNames,
    loader,
    NSigmasISysts,
    NSigmas,
    kTruthCut_IsSignal,
    kNoShift,
    true
  );


  // Run

  loader.Go();

  // merge trees

  tree_SelectedEvents_NSigmas->MergeTree( *tree_SelectedEvents );
  tree_Truth_NSigmas->MergeTree( *tree_Truth );


  // Output

  TFile *f_out = new TFile("MCTree.root", "RECREATE");
  tree_SelectedEvents_NSigmas->SaveToTClonesArrays(f_out);
  tree_Truth_NSigmas->SaveTo(f_out);

}
