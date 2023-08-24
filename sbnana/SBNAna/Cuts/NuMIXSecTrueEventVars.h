#pragma once

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/PrimaryUtils.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

namespace ana{

  std::vector<double> GetTrueVarVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    std::function<double(const caf::SRTrueInteractionProxy&) > trueth_var
  );

  std::vector<double> GetCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  );

  // Interaction
  const SpillMultiVar kCutTypeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetCutTypeVectorPerNu(sr, Is1muNp0piWithProtonPcut);
    }
  );
  const SpillMultiVar kTruePDGVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoPDG_True);
    }
  );
  const SpillMultiVar kTrueModeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoMode_True);
    }
  );
  const SpillMultiVar kTrueTargetVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::Target_True);
    }
  );
  const SpillMultiVar kTrueEVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoE_True);
    }
  );

  // Muon
  const SpillMultiVar kTrueMuonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::MuonP_True);
    }
  );
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::MuonNuCosineTheta_True);
    }
  );

  // Proton
  const SpillMultiVar kTrueProtonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonP_True);
    }
  );
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonNuCosineTheta_True);
    }
  );

  // Muon+Proton

  const SpillMultiVar kTrueCosThMuonProtonVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::CosThMuonProton_True);
    }
  );

  // TKI
  const SpillMultiVar kTruedeltaPTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPT_True);
    }
  );
  const SpillMultiVar kTruedeltaPTxVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPTx_True);
    }
  );
  const SpillMultiVar kTruedeltaPTyVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPTy_True);
    }
  );
  const SpillMultiVar kTruedeltaalphaTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaalphaT_True);
    }
  );
  const SpillMultiVar kTruedeltaphiTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaphiT_True);
    }
  );

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

}
