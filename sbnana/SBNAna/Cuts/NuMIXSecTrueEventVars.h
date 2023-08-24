#pragma once

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/PrimaryUtils.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

using namespace ana::PrimaryUtil;
using namespace caf;

namespace ana{

  std::vector<double> GetTrueVarVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const SRTrueInteractionProxy&)> isSignal,
    std::function<double(const SRTrueInteractionProxy&) > trueth_var
  );

  std::vector<double> GetCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const SRTrueInteractionProxy&)> isSignal
  );

  // Interaction
  const SpillMultiVar kCutTypeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetCutTypeVectorPerNu(sr, Is1muNp0piWithProtonPcut);
    }
  );
  const SpillMultiVar kTruePDGVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, NeutrinoPDG);
    }
  );
  const SpillMultiVar kTrueModeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, NeutrinoMode);
    }
  );
  const SpillMultiVar kTrueTargetVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, Target);
    }
  );
  const SpillMultiVar kTrueEVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, NeutrinoE);
    }
  );

  // Muon
  const SpillMultiVar kTrueMuonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, MuonP);
    }
  );
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, MuonNuCosineTheta);
    }
  );

  // Proton
  const SpillMultiVar kTrueProtonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ProtonP);
    }
  );
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ProtonNuCosineTheta);
    }
  );

  // Muon+Proton

  const SpillMultiVar kTrueCosThMuonProtonVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, CosThMuonProton);
    }
  );

  // TKI
  const SpillMultiVar kTruedeltaPTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, deltaPT);
    }
  );
  const SpillMultiVar kTruedeltaPTxVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, deltaPTx);
    }
  );
  const SpillMultiVar kTruedeltaPTyVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, deltaPTy);
    }
  );
  const SpillMultiVar kTruedeltaalphaTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, deltaalphaT);
    }
  );
  const SpillMultiVar kTruedeltaphiTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, deltaphiT);
    }
  );

  // Weight

  std::vector<double> GetNuMIPPFXWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const SRTrueInteractionProxy&)> isSignal
  );

  const SpillMultiVar kNuMIPPFXWeightVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetNuMIPPFXWeightVectorPerNu(sr, Is1muNp0piWithProtonPcut);
    }
  );

}
