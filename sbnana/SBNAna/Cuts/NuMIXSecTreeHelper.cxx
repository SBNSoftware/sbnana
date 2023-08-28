#include "sbnana/SBNAna/Cuts/NuMIXSecTreeHelper.h"

using namespace ana;
using namespace ana::PrimaryUtil;

namespace ana{

  std::vector<string> GetNuMITrueTreeLabels(){

    return {
      // CutType
      "CutType/i",
      // Interaction
      "TruePDG/i", "TrueMode/i", "TrueTarget/i",
      // Weight
      "FluxWeight",
      // Nu E
      "TrueE",
      // Muon
      "TrueMuonP",
      "TrueMuonPt",
      "TrueMuonCos",
      "TrueMuonLength",
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
  std::vector<SpillMultiVar> GetNuMITrueTreeVars(){

    return {
      // CutType
      kCutTypeVectorPerSignalNu,
      // Interaction
      kTruePDGVectorPerSignalNu, kTrueModeVectorPerSignalNu, kTrueTargetVectorPerSignalNu,
      // Weight
      kNuMIPPFXWeightVectorPerSignalNu,
      // Nu E
      kTrueEVectorPerSignalNu,
      // Muon
      kTrueMuonPVectorPerSignalNu,
      kTrueMuonPtVectorPerSignalNu,
      kTrueMuonNuCosineThetaVectorPerSignalNu,
      kTrueMuonLengthVectorPerSignalNu,
      // Proton
      kTrueProtonPVectorPerSignalNu,
      kTrueProtonPtVectorPerSignalNu,
      kTrueProtonNuCosineThetaVectorPerSignalNu,
      kTrueProtonLengthVectorPerSignalNu,
      // Muon+Proton
      kTrueCosThMuonProtonVectorPerSignalNu,
      // TKI
      kTruedeltaPTVectorPerSignalNu,
      kTruedeltaPTxVectorPerSignalNu,
      kTruedeltaPTyVectorPerSignalNu,
      kTruedeltaalphaTVectorPerSignalNu,
      kTruedeltaphiTVectorPerSignalNu,
    };

  }

  std::vector<std::string> GetNuMIRecoTreeLabels(){

    return {
      // CutType
      "CutType/i",
      // Intercation
      "TruePDG/i", "TrueMode/i", "TrueTarget/i", "TrueIsCC/i",
      // Number of primary particles
      "TrueNProton/i", "TrueNNeutron/i",
      "TrueNpip/i", "TrueNpim/i", "TrueNpi0/i",
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
    };

  }
  std::vector<Var> GetNuMIRecoTreeVars(){

    return {
      // CutType
      kNuMICutType,
      // Intercation
      kNuMITruePDG, kNuMITrueMode, kNuMITrueTarget, kNuMITrueIsCC,
      // Number of primary particles
      kNuMITrueNProton, kNuMITrueNNeutron,
      kNuMITrueNpip, kNuMITrueNpim, kNuMITrueNpi0,
      // Is Signal
      kNuMIIsSignal,
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
    };

  }

}
