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
      "TrueMuonCos",
      // Proton
      "TrueProtonP",
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
      kTrueMuonNuCosineThetaVectorPerSignalNu,
      // Proton
      kTrueProtonPVectorPerSignalNu,
      // Muon+Proton
      kTrueProtonNuCosineThetaVectorPerSignalNu,
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
      // Number of primary pions
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
      // Muon angle
      "RecoMuonCos", "TrueMuonCos",
      // Proton momentum
      "RecoProtonP", "TrueProtonP",
      // Proton angle
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
      // Number of primary pions
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
      // Muon angle
      kNuMIRecoCosThBeam, kNuMITrueCosThBeam,
      // Proton momentum
      kNuMIProtonCandidateRecoP, kNuMIProtonTrueP,
      // Proton angle
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
