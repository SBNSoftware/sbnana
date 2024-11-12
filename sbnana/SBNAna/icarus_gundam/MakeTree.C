#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
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
  vec_inputs.push_back("root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/archive/sam_managed_users/jskim/data/4/0/0/d/5a420bc4-9fd8-43cd-b7d0-4dbf5dc50127-flatcaf_SBNNuSyst_0.root");
  bool IsData = false;

  // SpectrumLoader
  SpectrumLoader loader(vec_inputs); //, kBeam, 10);

  // Define Tree variables
  // Example here are defined at sbnana/SBNAna/Vars/NuMIXSecVars.h
  std::vector<std::string> vec_labels_SelectedEvents = {
    "CutType/i",
    "IsSignal/i",
    "RecoMuonP", "TrueMuonP",
    "RecoMuonCos", "TrueMuonCos",
    "LeadingChargedPionCandidateLength",
  };
  std::vector<Var> vec_vars_SelectedEvents = {
    kNuMICutType,
    kNuMISliceSignalType,
    kNuMIMuonCandidateRecoP, kNuMIMuonTrueP,
    kNuMIRecoCosThVtx, kNuMITrueCosThVtx,
    kNuMILeadingChargedPionCandidateLength,
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

  std::vector<std::string>  GENIEMorphKnobNames = ICARUSGUNDAM::GetGENIEMorphKnobNames();
  for(unsigned int i=0; i<GENIEMorphKnobNames.size(); i++){
    NSigmasPsetNames.push_back( GENIEMorphKnobNames.at(i) );

    std::string psetname = SystProviderPrefix+"_multisigma_"+GENIEMorphKnobNames.at(i);
    NSigmasISysts.push_back( new SBNWeightMirrorSyst(psetname) );

    NSigmas.push_back( {-1, -0.5, 0, 0.5, 1} );
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

  // Run

  loader.Go();

  // merge trees

  tree_SelectedEvents_NSigmas->MergeTree( *tree_SelectedEvents );


  // Output

  TFile *f_out = new TFile("output.root", "RECREATE");
  tree_SelectedEvents_NSigmas->SaveToTClonesArrays(f_out);

}
