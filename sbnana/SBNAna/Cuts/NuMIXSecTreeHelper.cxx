#include "sbnana/SBNAna/Cuts/NuMIXSecTreeHelper.h"

using namespace ana;

namespace ana{

  std::vector<string> GetNuMITrueTreeLabels(){

    return {
      // Interaction
      "TruePDG/i", "TrueMode/i", "TrueTarget/i", "TrueIsCC/i",
      "TrueIsFHC/i",
      // Weight
      "FluxWeight",
      // Nu E
      "TrueE",
      // Muon
      "TrueMuonP",
      "TrueMuonPt",
      "TrueMuonCos",
      "TrueMuonCosBeam",
      "TrueMuonLength",
      "TrueMuonContained/i",
      // Proton
      "TrueProtonP",
      "TrueProtonPt",
      "TrueProtonCos",
      "TrueProtonLength",
      // Muon+Proton
      "TrueMuonProtonCos",
      // TKI
      "TruedeltaPT",
      "TruedeltaPTx",
      "TruedeltaPTy",
      "TruedeltaalphaT",
      "TruedeltaphiT",
    };

  }
  std::vector<TruthVar> GetNuMITrueTreeVars(){

    return {
      // Interaction
      kTruth_NeutrinoPDG, kTruth_NeutrinoMode, kTruth_Target, kTruth_IsCC,
      kTruth_IsFHC,
      // Weight
      kGetTruthNuMIFluxWeight,
      // Nu E
      kTruth_NeutrinoE,
      // Muon
      kTruth_MuonP,
      kTruth_MuonPt,
      kTruth_MuonNuCosineTheta,
      kTruth_MuonCosThBeam,
      kTruth_MuonLength,
      kTruth_MuonContained,
      // Proton
      kTruth_ProtonP,
      kTruth_ProtonPt,
      kTruth_ProtonNuCosineTheta,
      kTruth_ProtonLength,
      // Muon+Proton
      kTruth_CosThMuonProton,
      // TKI
      kTruth_deltaPT,
      kTruth_deltaPTx,
      kTruth_deltaPTy,
      kTruth_deltaalphaT,
      kTruth_deltaphiT,
    };

  }

  std::vector<std::string> GetNuMIRecoTreeLabels(){

    return {
      // CutType
      "CutType/i",
      // Intercation
      "TruePDG/i", "TrueMode/i", "TrueTarget/i", "TrueIsCC/i",
      "TrueIsFHC/i",
      // Number of primary particles
      "TrueNProton/i", "TrueNNeutron/i",
      "TrueNpip/i", "TrueNpim/i", "TrueNpi0/i",
      "TrueNpipAll/i", "TrueNpimAll/i", "TrueNpi0All/i",
      // True muon
      "TrueMuonContained/i",
      // Is Signal
      "IsSignal/i",
      // Weight
      "FluxWeight",
      // NuE
      "TrueE",
      // RecoMuon track notfound/contained/exiting
      "MuonTrackType/i",
      // Muon momentum
      "RecoMuonP", "TrueMuonP",
      "RecoMuonPt", "TrueMuonPt",
      // Muon length
      "RecoMuonLength", "TrueMuonLength",
      // Muon angle
      "RecoMuonCos", "TrueMuonCos",
      "RecoMuonCosBeam", "TrueMuonCosBeam",
      // Proton momentum
      "RecoProtonP", "TrueProtonP",
      "RecoProtonPt", "TrueProtonPt",
      // Proton length
      "RecoProtonLength", "TrueProtonLength",
      // Proton angle
      "RecoProtonCos", "TrueProtonCos",
      // Muon,Proton angle
      "RecoMuonProtonCos", "TrueMuonProtonCos",
      // TKI
      "RecodeltaPT", "TruedeltaPT",
      "RecodeltaPTx", "TruedeltaPTx",
      "RecodeltaPTy", "TruedeltaPTy",
      "RecodeltaalphaT", "TruedeltaalphaT",
      "RecodeltaphiT", "TruedeltaphiT",
      // SideBand
      // pi+-
      // - True,
      "TrueChargedPionKE",
      // - Reco
      //   - Track
      "LeadingChargedPionCandidateLength",
      "LeadingChargedPionCandidateNDaughter/i",
      "LeadingChargedPionCandidateMatchedPDG/i",
      "LeadingChargedPionCandidateNCollectionHit/i",
      "LeadingChargedPionCandidateMIPChi2",
      // pi0
      "LeadingPhotonCandidateE",
      "SecondaryPhotonCandidateE",
      "PhotonCandidatesOpeningAngle",
    };

  }
  std::vector<Var> GetNuMIRecoTreeVars(){

    return {
      // CutType
      kNuMICutType,
      // Intercation
      kNuMITruePDG, kNuMITrueMode, kNuMITrueTarget, kNuMITrueIsCC,
      kNuMIIsFHC,
      // Number of primary particles
      kNuMITrueNProton, kNuMITrueNNeutron,
      kNuMITrueNpip, kNuMITrueNpim, kNuMITrueNpi0,
      kNuMITrueNpip_All, kNuMITrueNpim_All, kNuMITrueNpi0_All,
      // True muon
      kNuMITrueMuonContained,
      // Is Signal
      kNuMISliceSignalType,
      // Weight
      kGetNuMIFluxWeight,
      // NuE
      kNuMITrueNuE,
      // RecoMuon track notfound/contained/exiting
      kNuMIRecoMuonContained,
      // Muon momentum
      kNuMIMuonCandidateRecoP, kNuMIMuonTrueP,
      kNuMIRecoMuonPt, kNuMITrueMuonPt,
      // Muon length
      kNuMIRecoMuonLength, kNuMITrueMuonLength,
      // Muon angle
      kNuMIRecoCosThVtx, kNuMITrueCosThVtx,
      kNuMIRecoCosThBeam, kNuMITrueCosThBeam,
      // Proton momentum
      kNuMIProtonCandidateRecoP, kNuMIProtonTrueP,
      kNuMIRecoProtonPt, kNuMITrueProtonPt,
      // Proton legnth
      kNuMIRecoProtonLength, kNuMITrueProtonLength,
      // Proton angle
      kNuMIProtonRecoCosThVtx, kNuMIProtonTrueCosThVtx,
      // Muon,Proton angle
      kNuMIRecoCosThMuP, kNuMITrueCosThMuP,
      // TKI
      kNuMIRecodeltaPT, kNuMITruedeltaPT,
      kNuMIRecodeltaPTx, kNuMITruedeltaPTx,
      kNuMIRecodeltaPTy, kNuMITruedeltaPTy,
      kNuMIRecodeltaalphaT, kNuMITruedeltaalphaT,
      kNuMIRecodeltaphiT, kNuMITruedeltaphiT,
      // SideBand
      // pi+-
      // - True
      kNuMITrueChargedPionKE,
      // - Reco
      //   - Track
      kNuMILeadingChargedPionCandidateLength,
      kNuMILeadingChargedPionCandidateNDaughter,
      kNuMILeadingChargedPionCandidateMatchedPDG,
      kNuMILeadingChargedPionCandidateNCollectionHit,
      kNuMILeadingChargedPionCandidateMIPChi2,
      // pi0
      kNuMILeadingPhotonCandidateE,
      kNuMISecondaryPhotonCandidateE,
      kNuMIPhotonCandidatesOpeningAngle,
    };

  }

}
