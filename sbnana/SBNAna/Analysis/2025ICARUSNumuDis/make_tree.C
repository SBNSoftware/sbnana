#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/SBNOnOffSysts.h"
#include "sbnana/CAFAna/Systs/UniverseOracle.h"
#include "sbnana/CAFAna/Systs/IcarusRun2DetectorSysts.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202401.h"
#include "sbnana/SBNAna/Cuts/ICARUSDataQualityCuts.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/SBNAna/Vars/Vars.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TDirectory.h"
#include "TFile.h"

#include <string>
#include <utility>
#include <vector>

using namespace ana;

std::vector<std::string> bound_at_neg2{
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NormCCMEC",
        "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
};

double offbeam_livetime = 0;
const SpillVar kOffbeamLivetime([](const caf::SRSpillProxy *sr) {
  if(icarus::kGoodRunsRun2(sr))
    offbeam_livetime += sr->hdr.noffbeambnb;
  return 1;
});

void make_tree(std::string outname = "/exp/icarus/data/users/jlarkin/icarus_numudis_sbruce.root")
{
  SpectrumLoader mc("/pnfs/icarus/persistent/users/jlarkin/cafs/icarus_numudis_cv_*.flat.caf.root");
  SpectrumLoader offbeam("/pnfs/icarus/persistent/users/jlarkin/cafs/run2_offbeam_*.flat.caf.root");
  
  const Var kTrueE = SIMPLEVAR(truth.E);
  const Var kTrueL = SIMPLEVAR(truth.baseline);
  const Var kTruePDG = SIMPLEVAR(truth.pdg);
  const Var kTrueCC = SIMPLEVAR(truth.iscc);
  const Var kSlcVX = SIMPLEVAR(vertex.x);
  const Var kSlcVY = SIMPLEVAR(vertex.y);
  const Var kSlcVZ = SIMPLEVAR(vertex.z);

  const Cut kTrueNu = SIMPLEVAR(truth.index) >= 0;

  const SpillCut kSpillSelection = kIcarus202401CRTPMTVeto;
  const Cut kSliceSelection = kIcarus202401Contained1muNp;

  std::vector<std::string> nu_branch_names = {
    "trueE", "trueL", "truePDG", "CC", "trueLeadingPMom", "trueLeadingPCosTheta", "trueMuMom", "trueMuCosTheta",
    "recoE", "reco_muon_p", "reco_muon_angle","deltaPt","slcVtxX", "slcVtxY", "slcVtxZ","endMuonX", "endMuonY", "endMuonZ", "deltaZ", "reco_muon_phi",
    "leading_proton_len", "leading_proton_p", "muon_length", "muon_dirx", "muon_diry", "muon_dirz", "nproton",
    "leading_proton_dirx", "leading_proton_diry", "leading_proton_dirz", "angle_mup", "chi2_mu", "chi2_p",
    "leading_proton_endx", "leding_proton_endy", "leading_proton_endz",
  };

  std::vector<Var> nu_vars = {
    kTrueE, kTrueL, kTruePDG, kTrueCC,kICARUS202401TrueLeadingProtonMomentum, 
    kICARUS202401TrueLeadingProtonCosTheta, kIcarus202401TrueMuonP, kIcarus202401TrueMuonCosTheta,
    kIcarus202401RecoENu, kIcarus202401RecoMuonP, kIcarus202401RecoMuonCosTheta, kIcarus202401RecoTransP, 
    kSlcVX, kSlcVY, kSlcVZ, kIcarus202401RecoMuonEndX, kIcarus202401RecoMuonEndY, 
    kIcarus202401RecoMuonEndZ, kIcarus202401BaryFMDeltaZ, kIcarus202401RecoMuonPhi, 
    kIcarus202401RecoLeadingProtonLen, kIcarus202401RecoLeadingProtonP, kIcarus202401RecoMuonLen, 
    kIcarus202401RecoMuonDirX, kIcarus202401RecoMuonDirY, kIcarus202401RecoMuonDirZ, 
    kIcarus202401NumProtons, kIcarus202401RecoLeadingProtonDirX, kIcarus202401RecoLeadingProtonDirY, 
    kIcarus202401RecoLeadingProtonDirZ, kIcarus202401RecoMuLeadPOpeningAngle, kIcarus202401MuonChi2Mu, 
    kIcarus202401LeadingProtonChi2Proton, kIcarus202401RecoLeadingProtonEndX, 
    kIcarus202401RecoLeadingProtonEndY, kIcarus202401RecoLeadingProtonEndZ, 
  };

  std::vector<std::string> branch_names = {
    "recoE", "reco_muon_p", "reco_muon_angle","deltaPt","slcVtxX", "slcVtxY", "slcVtxZ","endMuonX", "endMuonY", "endMuonZ", "deltaZ", "reco_muon_phi",
    "leading_proton_len", "leading_proton_p", "muon_length", "muon_dirx", "muon_diry", "muon_dirz", "nproton",
    "leading_proton_dirx", "leading_proton_diry", "leading_proton_dirz", "angle_mup", "chi2_mu", "chi2_p",
    "leading_proton_endx", "leding_proton_endy", "leading_proton_endz",
  };

  std::vector<Var> vars = {
    kIcarus202401RecoENu, kIcarus202401RecoMuonP, kIcarus202401RecoMuonCosTheta, kIcarus202401RecoTransP, 
    kSlcVX, kSlcVY, kSlcVZ, kIcarus202401RecoMuonEndX, kIcarus202401RecoMuonEndY, 
    kIcarus202401RecoMuonEndZ, kIcarus202401BaryFMDeltaZ, kIcarus202401RecoMuonPhi, 
    kIcarus202401RecoLeadingProtonLen, kIcarus202401RecoLeadingProtonP, kIcarus202401RecoMuonLen, 
    kIcarus202401RecoMuonDirX, kIcarus202401RecoMuonDirY, kIcarus202401RecoMuonDirZ, 
    kIcarus202401NumProtons, kIcarus202401RecoLeadingProtonDirX, kIcarus202401RecoLeadingProtonDirY, 
    kIcarus202401RecoLeadingProtonDirZ, kIcarus202401RecoMuLeadPOpeningAngle, kIcarus202401MuonChi2Mu,
    kIcarus202401LeadingProtonChi2Proton, kIcarus202401RecoLeadingProtonEndX, 
    kIcarus202401RecoLeadingProtonEndY, kIcarus202401RecoLeadingProtonEndZ, 
  };

  Tree nutree("selectedNu", nu_branch_names, mc, nu_vars, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);
  Tree costree("selectedCos", branch_names, mc, vars, kSpillSelection, kSliceSelection && !kTrueNu, kNoShift, true, true);
  Tree offbeamtree("selectedOffbeam", branch_names, offbeam, vars, kSpillSelection && icarus::kGoodRunsRun2, kSliceSelection, kNoShift, true, true);
  Spectrum dummy_spec("", Binning::Simple(2,0,2), offbeam, kOffbeamLivetime, kNoSpillCut); 

  std::vector<std::string> genie_names = GetSBNGenieWeightNames();
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b1");
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b2");
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b3");
  genie_names.push_back("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b4");
  genie_names.push_back("GENIEReWeight_SBNNuSyst_GENIE_multisigma_FracPN_CCMEC");
  SBNWeightSyst b1("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b1");
  SBNWeightSyst b2("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b2");
  SBNWeightSyst b3("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b3");
  SBNWeightSyst b4("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b4");
  SBNWeightSyst fracPN("GENIEReWeight_SBNNuSyst_GENIE_multisigma_FracPN_CCMEC");
  std::vector<const ISyst*> genie_systs = GetSBNGenieWeightSysts();
  genie_systs.push_back(&b1);
  genie_systs.push_back(&b2);
  genie_systs.push_back(&b3);
  genie_systs.push_back(&b4);
  genie_systs.push_back(&fracPN);
  std::vector<std::string> onoff_names = GetSBNOnOffNames();
  std::vector<const ISyst*> onoff_systs = GetSBNOnOffSysts();
    onoff_names.push_back("GENIEReWeight_SBNNuSyst_GENIE_multisigma_XSecShape_CCMEC");
  onoff_names.push_back("FSIReweight_SBNNuSyst_FSI_hNReweight_multisigma_FSIReweight");
  onoff_names.push_back("FSIReweight_SBNNuSyst_FSI_INCLReweight_multisigma_FSIReweight");
  onoff_names.push_back("FSIReweight_SBNNuSyst_FSI_G4BCReweight_multisigma_FSIReweight");
  SBNOnOffSyst mecshape("GENIEReWeight_SBNNuSyst_GENIE_multisigma_XSecShape_CCMEC");
  SBNOnOffSyst hN_FSI("FSIReweight_SBNNuSyst_FSI_hNReweight_multisigma_FSIReweight");
  SBNOnOffSyst INCL_FSI("FSIReweight_SBNNuSyst_FSI_INCLReweight_multisigma_FSIReweight");
  SBNOnOffSyst G4BC_FSI("FSIReweight_SBNNuSyst_FSI_G4BCReweight_multisigma_FSIReweight");
  onoff_systs.push_back(&mecshape);
  onoff_systs.push_back(&hN_FSI);
  onoff_systs.push_back(&INCL_FSI);
  onoff_systs.push_back(&G4BC_FSI);

  std::vector<std::string> nsigma_names = genie_names;
  std::vector<const ISyst*> nsigma_systs = genie_systs;
  for(const std::string &name: onoff_names) nsigma_names.push_back(name);
  for(const ISyst *s: onoff_systs) nsigma_systs.push_back(s);
  std::vector<std::pair<int, int>> min_max;
  for(size_t i = 0; i < genie_names.size(); ++i) min_max.emplace_back(-3, 3);
  for(size_t i = 0; i < onoff_names.size(); ++i) min_max.emplace_back(0, 1);
  for(size_t i = 0; i < genie_names.size(); ++i) {
    if(std::find(bound_at_neg2.begin(), bound_at_neg2.end(), genie_names[i]) != bound_at_neg2.end())
      min_max[i] = std::make_pair(-2, 3);
  }

  NSigmasTree nsigtree("multisigmaTree", nsigma_names, mc, nsigma_systs, min_max, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);

  const std::vector<std::string> flux_names{ "expskin_Flux", "horncurrent_Flux", "nucleoninexsec_Flux", "nucleonqexsec_Flux", "nucleontotxsec_Flux", "pioninexsec_Flux", "pionqexsec_Flux", "piontotxsec_Flux", "piplus_Flux", "piminus_Flux", "kplus_Flux", "kminus_Flux", "kzero_Flux" };
  const std::vector<std::string> xsec_multisim_names{
    "GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_RPA_CCQE",
    "GENIEReWeight_SBN_v1_multisim_CoulombCCQE",
    "GENIEReWeight_SBN_v1_multisim_NormCCMEC",
    "GENIEReWeight_SBN_v1_multisim_NormNCMEC",
    "GENIEReWeight_SBN_v1_multisim_NCELVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi",
    "GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi",
    "GENIEReWeight_SBN_v1_multisim_COHVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse",
    "GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse",
  };
  std::vector<std::string> multisim_names;
  std::vector<std::vector<Var>> univsKnobs;
  std::vector<unsigned int> nuniverses;
  for(const auto& name: flux_names) {
    multisim_names.push_back(name);
    size_t nuniv = 1000;
    nuniverses.push_back(nuniv);
    univsKnobs.emplace_back();
    for(size_t i = 0; i < nuniv; ++i) 
      univsKnobs.back().push_back(GetUniverseWeight(name, i));
  }

  for(const auto& name: xsec_multisim_names) {
    multisim_names.push_back(name);
    size_t nuniv = 100;
    nuniverses.push_back(nuniv);
    univsKnobs.emplace_back();
    for(size_t i = 0; i < nuniv; ++i) 
      univsKnobs.back().push_back(GetUniverseWeight(name, i));
  }
  multisim_names.push_back("GENIEReWeight_SBNNuSyst_GENIE_multisigma_DecayAngMECVariationResponse");
  nuniverses.push_back(25);
  univsKnobs.emplace_back();
  for(size_t i = 0; i < 25; ++i)
    univsKnobs.back().push_back(GetUniverseWeight("GENIEReWeight_SBNNuSyst_GENIE_multisigma_DecayAngMECVariationResponse", i));
  multisim_names.push_back("GENIEReWeight_SBNNuSyst_LQCDZExpFit_multisim_ZExpAVariationResponse");
  nuniverses.push_back(2);
  univsKnobs.emplace_back();
  univsKnobs.back().push_back(GetUniverseWeight("GENIEReWeight_SBNNuSyst_LQCDZExpFit_multisim_ZExpAVariationResponse", 0));
  univsKnobs.back().push_back(GetUniverseWeight("GENIEReWeight_SBNNuSyst_LQCDZExpFit_multisim_ZExpAVariationResponse", 1));
  
  std::string name1 = "deltapT";
  std::string name2 = "recoE";
  std::string name3 = "L_p";
  std::string name4 = "L_mu";

  std::vector<const ISyst*> detsysts  = GetIcarusRun2DetectorSysts(name1, kIcarus202401RecoTransP);
  std::vector<const ISyst*> detsysts2 = GetIcarusRun2DetectorSysts(name2, kIcarus202401RecoENu);
  std::vector<const ISyst*> detsysts3 = GetIcarusRun2DetectorSysts(name3, kIcarus202401RecoLeadingProtonLen);
  std::vector<const ISyst*> detsysts4 = GetIcarusRun2DetectorSysts(name4, kIcarus202401RecoMuonLen);
 
  for(const auto& syst: detsysts2) detsysts.push_back(syst);
  for(const auto& syst: detsysts3) detsysts.push_back(syst);
  for(const auto& syst: detsysts4) detsysts.push_back(syst);

  std::vector<std::string> detsyst_names;
  std::vector<std::pair<int,int>> min_max_det;
  for(size_t i = 0; i < detsysts.size(); ++i) {
    min_max_det.emplace_back(-3,3);
    const ISyst* detsyst = detsysts[i];
    detsyst_names.push_back(detsyst->ShortName()+"_multisigma");
  }

  NSigmasTree detsysttree("detsystTree", detsyst_names, mc, detsysts, min_max_det, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);

  NUniversesTree nunivtree("multisimTree", multisim_names, mc, univsKnobs, nuniverses, kSpillSelection, kSliceSelection && kTrueNu, kNoShift, true, true);

  mc.Go();
  offbeam.Go();

  offbeamtree.OverrideLivetime(offbeam_livetime);

  //nsigtree.MergeTree(nutree);
  //nunivtree.MergeTree(nsigtree);

  TFile fout(outname.c_str(), "RECREATE");
  TDirectory* dir = fout.mkdir("events");
  nutree.SaveTo(dir); 
  costree.SaveTo(dir);
  nsigtree.SaveTo(dir);
  nunivtree.SaveTo(dir);
  detsysttree.SaveTo(dir);
  TDirectory *offdir = fout.mkdir("offbeam");
  offbeamtree.SaveTo(offdir);
}
