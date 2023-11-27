#pragma once

#include "sbnana/CAFAna/Systs/UniverseOracle.h"

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/PrimaryUtils.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

namespace ana{

  // Return a vector of TrueVar for each true neutrino that passes isSignal
  std::vector<double> GetTrueVarVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    std::function<double(const caf::SRTrueInteractionProxy&) > TrueVar
  );
  // Use GetTrueVarVectorPerNu and create a SpillMultiVar
  SpillMultiVar GetTrueSpillMultiVarPerSignalNu(
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    std::function<double(const caf::SRTrueInteractionProxy&) > TrueVar
  );

  // For each true signal neutrino,
  // return 1 if it is matched to a reco slice that pass signal selection
  // return 0 if it is NOT matched dto any reco slice that pass signal selection
  std::vector<double> GetCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  );
  // Use GetCutTypeVectorPerNu and create a SpillMultiVar
  const SpillMultiVar kCutTypeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetCutTypeVectorPerNu(sr, Is1muNp0piWithPhaseSpaceCut);
    }
  );

  // Similar but checks against the SpillCut
  std::vector<double> GetSpillCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  );
  // Use GetSpillCutTypeVectorPerNu and create a SpillMultiVar
  const SpillMultiVar kSpillCutTypeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetSpillCutTypeVectorPerNu(sr, Is1muNp0piWithPhaseSpaceCut);
    }
  );
  const SpillMultiVar kSpillCutTypeVectorPerFVNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetSpillCutTypeVectorPerNu(sr, IsNuInFV);
    }
  );
  const SpillMultiVar kSpillCutTypeVectorPerAVNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetSpillCutTypeVectorPerNu(sr, IsNuInAV);
    }
  );

  // For each true signal neutrino,
  // return 1 if it is matched to a reco slice that pass signal selection (WITHOUT THE SHOWER CUT)
  // return 0 if it is NOT matched dto any reco slice that pass signal selection
  std::vector<double> GetCutTypeWithoutShowerCutVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  );
  // Use GetCutTypeVectorPerNu and create a SpillMultiVar
  const SpillMultiVar kCutTypeWithoutShowerCutVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetCutTypeWithoutShowerCutVectorPerNu(sr, Is1muNp0piWithPhaseSpaceCut);
    }
  );

  // Interaction
  const SpillMultiVar kTruePDGVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::NeutrinoPDG_True);
  const SpillMultiVar kTrueModeVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::NeutrinoMode_True);
  const SpillMultiVar kTrueIsCCVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::IsCC_True);
  const SpillMultiVar kTrueTargetVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::Target_True);
  const SpillMultiVar kTrueEVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::NeutrinoE_True);
  const SpillMultiVar kTrueQ2VectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::Q2_True);
  const SpillMultiVar kTrueq0VectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::q0_True);
  const SpillMultiVar kTrueq3VectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::q3_True);
  const SpillMultiVar kTruewVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::w_True);

  // Muon
  const SpillMultiVar kTrueMuonPVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::MuonP_True);
  const SpillMultiVar kTrueMuonPtVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::MuonPt_True);
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::MuonNuCosineTheta_True);
  const SpillMultiVar kTrueMuonCosThBeamVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::MuonCosThBeam_True);
  const SpillMultiVar kTrueMuonLengthVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::MuonLength_True);
  const SpillMultiVar kTrueMuonContainedVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::MuonContained_True);

  // Proton
  const SpillMultiVar kTrueProtonPVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::ProtonP_True);
  const SpillMultiVar kTrueProtonPtVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::ProtonPt_True);
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::ProtonNuCosineTheta_True);
  const SpillMultiVar kTrueProtonLengthVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::ProtonLength_True);

  // Muon+Proton
  const SpillMultiVar kTrueCosThMuonProtonVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::CosThMuonProton_True);

  // TKI
  const SpillMultiVar kTruedeltaPTVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::deltaPT_True);
  const SpillMultiVar kTruedeltaPTxVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::deltaPTx_True);
  const SpillMultiVar kTruedeltaPTyVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::deltaPTy_True);
  const SpillMultiVar kTruedeltaalphaTVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::deltaalphaT_True);
  const SpillMultiVar kTruedeltaphiTVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithPhaseSpaceCut, PrimaryUtil::deltaphiT_True);

  // Weight
  std::vector<double> GetNuMIPPFXWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  );
  const SpillMultiVar kNuMIPPFXWeightVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetNuMIPPFXWeightVectorPerNu(sr, Is1muNp0piWithPhaseSpaceCut);
    }
  );
  const SpillMultiVar kNuMIPPFXWeightVectorPerFVNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetNuMIPPFXWeightVectorPerNu(sr, IsNuInFV);
    }
  );
  const SpillMultiVar kNuMIPPFXWeightVectorPerAVNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetNuMIPPFXWeightVectorPerNu(sr, IsNuInAV);
    }
  );

  std::vector<double> GetSigmaWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    const std::string& psetName,
    double shift // in sigma
  );
  SpillMultiVar GetSigmaWeightSpillMultiVarPerSignalNu(
    const std::string& psetName,
    double shift // in sigma
  );

  struct Univs
  {
    int i0, i1;
    double w0, w1;
  };

  std::vector<double> GetUniverseWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    const std::string& psetName,
    int univIdx
  );
  SpillMultiVar GetUniverseWeightSpillMultiVarPerSignalNu(
    const std::string& psetName,
    int univIdx
  );

  /// And now many of the same as for 0pi signal, but for any neutrino in the fiducial volume (IsNuInFV)
  const SpillMultiVar kTruePDGVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::NeutrinoPDG_True);
  const SpillMultiVar kTrueModeVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::NeutrinoMode_True);
  const SpillMultiVar kTrueIsCCVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::IsCC_True);
  const SpillMultiVar kTrueTargetVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::Target_True);
  const SpillMultiVar kTrueEVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::NeutrinoE_True);
  const SpillMultiVar kTrueQ2VectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::Q2_True);
  const SpillMultiVar kTrueq0VectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::q0_True);
  const SpillMultiVar kTrueq3VectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::q3_True);
  const SpillMultiVar kTruewVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::w_True);
  const SpillMultiVar kTrueMuonPVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::MuonP_True);
  const SpillMultiVar kTrueMuonPtVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::MuonPt_True);
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::MuonNuCosineTheta_True);
  const SpillMultiVar kTrueMuonCosThBeamVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::MuonCosThBeam_True);
  const SpillMultiVar kTrueMuonLengthVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::MuonLength_True);
  const SpillMultiVar kTrueMuonContainedVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::MuonContained_True);
  const SpillMultiVar kTrueProtonPVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::ProtonP_True);
  const SpillMultiVar kTrueProtonPtVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::ProtonPt_True);
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::ProtonNuCosineTheta_True);
  const SpillMultiVar kTrueProtonLengthVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::ProtonLength_True);
  const SpillMultiVar kTrueCosThMuonProtonVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::CosThMuonProton_True);
  const SpillMultiVar kTruedeltaPTVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::deltaPT_True);
  const SpillMultiVar kTruedeltaPTxVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::deltaPTx_True);
  const SpillMultiVar kTruedeltaPTyVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::deltaPTy_True);
  const SpillMultiVar kTruedeltaalphaTVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::deltaalphaT_True);
  const SpillMultiVar kTruedeltaphiTVectorPerFVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInFV, PrimaryUtil::deltaphiT_True);

  // And active volume (AV)
  const SpillMultiVar kTruePDGVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::NeutrinoPDG_True);
  const SpillMultiVar kTrueModeVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::NeutrinoMode_True);
  const SpillMultiVar kTrueIsCCVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::IsCC_True);
  const SpillMultiVar kTrueTargetVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::Target_True);
  const SpillMultiVar kTrueIsInFVPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::IsInFV_True);
  const SpillMultiVar kTrueEVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::NeutrinoE_True);
  const SpillMultiVar kTrueQ2VectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::Q2_True);
  const SpillMultiVar kTrueq0VectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::q0_True);
  const SpillMultiVar kTrueq3VectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::q3_True);
  const SpillMultiVar kTruewVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::w_True);
  const SpillMultiVar kTrueMuonPVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::MuonP_True);
  const SpillMultiVar kTrueMuonPtVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::MuonPt_True);
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::MuonNuCosineTheta_True);
  const SpillMultiVar kTrueMuonCosThBeamVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::MuonCosThBeam_True);
  const SpillMultiVar kTrueMuonLengthVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::MuonLength_True);
  const SpillMultiVar kTrueMuonContainedVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::MuonContained_True);
  const SpillMultiVar kTrueProtonPVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::ProtonP_True);
  const SpillMultiVar kTrueProtonPtVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::ProtonPt_True);
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::ProtonNuCosineTheta_True);
  const SpillMultiVar kTrueProtonLengthVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::ProtonLength_True);
  const SpillMultiVar kTrueCosThMuonProtonVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::CosThMuonProton_True);
  const SpillMultiVar kTruedeltaPTVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::deltaPT_True);
  const SpillMultiVar kTruedeltaPTxVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::deltaPTx_True);
  const SpillMultiVar kTruedeltaPTyVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::deltaPTy_True);
  const SpillMultiVar kTruedeltaalphaTVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::deltaalphaT_True);
  const SpillMultiVar kTruedeltaphiTVectorPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::deltaphiT_True);
  const SpillMultiVar kTrueNuBaselinePerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::nuBaseline_True);
  const SpillMultiVar kTrueNuParentDkXPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::nuParentDkX_True);
  const SpillMultiVar kTrueNuParentDkYPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::nuParentDkY_True);
  const SpillMultiVar kTrueNuParentDkZPerAVNu = GetTrueSpillMultiVarPerSignalNu(IsNuInAV, PrimaryUtil::nuParentDkZ_True);
}
