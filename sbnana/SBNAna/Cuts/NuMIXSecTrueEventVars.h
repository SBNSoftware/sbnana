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
      return GetCutTypeVectorPerNu(sr, Is1muNp0piWithProtonPcut);
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
      return GetSpillCutTypeVectorPerNu(sr, Is1muNp0piWithProtonPcut);
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
      return GetCutTypeWithoutShowerCutVectorPerNu(sr, Is1muNp0piWithProtonPcut);
    }
  );

  // Interaction
  const SpillMultiVar kTruePDGVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoPDG_True);
  const SpillMultiVar kTrueModeVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoMode_True);
  const SpillMultiVar kTrueIsCCVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::IsCC_True);
  const SpillMultiVar kTrueTargetVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::Target_True);
  const SpillMultiVar kTrueEVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoE_True);
  const SpillMultiVar kTrueQ2VectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::Q2_True);
  const SpillMultiVar kTrueq0VectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::q0_True);
  const SpillMultiVar kTrueq3VectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::q3_True);
  const SpillMultiVar kTruewVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::w_True);

  // Muon
  const SpillMultiVar kTrueMuonPVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::MuonP_True);
  const SpillMultiVar kTrueMuonPtVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::MuonPt_True);
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::MuonNuCosineTheta_True);
  const SpillMultiVar kTrueMuonCosThBeamVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::MuonCosThBeam_True);
  const SpillMultiVar kTrueMuonLengthVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::MuonLength_True);

  // Proton
  const SpillMultiVar kTrueProtonPVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonP_True);
  const SpillMultiVar kTrueProtonPtVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonPt_True);
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonNuCosineTheta_True);
  const SpillMultiVar kTrueProtonLengthVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonLength_True);

  // Muon+Proton
  const SpillMultiVar kTrueCosThMuonProtonVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::CosThMuonProton_True);

  // TKI
  const SpillMultiVar kTruedeltaPTVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPT_True);
  const SpillMultiVar kTruedeltaPTxVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPTx_True);
  const SpillMultiVar kTruedeltaPTyVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPTy_True);
  const SpillMultiVar kTruedeltaalphaTVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::deltaalphaT_True);
  const SpillMultiVar kTruedeltaphiTVectorPerSignalNu = GetTrueSpillMultiVarPerSignalNu(Is1muNp0piWithProtonPcut, PrimaryUtil::deltaphiT_True);

  // Weight
  std::vector<double> GetNuMIPPFXWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  );
  const SpillMultiVar kNuMIPPFXWeightVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetNuMIPPFXWeightVectorPerNu(sr, Is1muNp0piWithProtonPcut);
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

}
