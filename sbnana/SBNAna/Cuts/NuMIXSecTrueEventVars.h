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
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoPDG);
    }
  );
  const SpillMultiVar kTrueModeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoMode);
    }
  );
  const SpillMultiVar kTrueTargetVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::Target);
    }
  );
  const SpillMultiVar kTrueEVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::NeutrinoE);
    }
  );

  // Muon
  const SpillMultiVar kTrueMuonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::MuonP);
    }
  );
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::MuonNuCosineTheta);
    }
  );

  // Proton
  const SpillMultiVar kTrueProtonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonP);
    }
  );
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::ProtonNuCosineTheta);
    }
  );

  // Muon+Proton

  const SpillMultiVar kTrueCosThMuonProtonVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::CosThMuonProton);
    }
  );

  // TKI
  const SpillMultiVar kTruedeltaPTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPT);
    }
  );
  const SpillMultiVar kTruedeltaPTxVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPTx);
    }
  );
  const SpillMultiVar kTruedeltaPTyVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaPTy);
    }
  );
  const SpillMultiVar kTruedeltaalphaTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaalphaT);
    }
  );
  const SpillMultiVar kTruedeltaphiTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, PrimaryUtil::deltaphiT);
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
