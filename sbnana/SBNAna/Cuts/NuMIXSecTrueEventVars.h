#pragma once

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/PrimaryUtils.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

using namespace ana;
using namespace ana::PrimaryUtil;

namespace ana{

  std::vector<double> GetTrueVarVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const TrueInteraction&)> isSignal,
    std::function<double(const TrueInteraction&) > trueth_var
  );

  std::vector<double> GetCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const TrueInteraction&)> isSignal
  );

  // Interaction
  const SpillMultiVar kCutTypeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetCutTypeVectorPerNu(sr, Is1muNp0piWithProtonPcut);
    }
  );
  const SpillMultiVar kTruePDGVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::NeutrinoPDG);
    }
  );
  const SpillMultiVar kTrueModeVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::NeutrinoMode);
    }
  );
  const SpillMultiVar kTrueTargetVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::Target);
    }
  );
  const SpillMultiVar kTrueEVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::NeutrinoE);
    }
  );

  // Muon
  const SpillMultiVar kTrueMuonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::MuonP);
    }
  );
  const SpillMultiVar kTrueMuonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::MuonNuCosineTheta);
    }
  );

  // Proton
  const SpillMultiVar kTrueProtonPVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::ProtonP);
    }
  );
  const SpillMultiVar kTrueProtonNuCosineThetaVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::ProtonNuCosineTheta);
    }
  );

  // Muon+Proton

  const SpillMultiVar kTrueCosThMuonProtonVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::CosThMuonProton);
    }
  );

  // TKI
  const SpillMultiVar kTruedeltaPTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::deltaPT);
    }
  );
  const SpillMultiVar kTruedeltaPTxVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::deltaPTx);
    }
  );
  const SpillMultiVar kTruedeltaPTyVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::deltaPTy);
    }
  );
  const SpillMultiVar kTruedeltaalphaTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::deltaalphaT);
    }
  );
  const SpillMultiVar kTruedeltaphiTVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetTrueVarVectorPerNu(sr, Is1muNp0piWithProtonPcut, ana::PrimaryUtil::deltaphiT);
    }
  );

  // Weight

  std::vector<double> GetNuMIPPFXWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const TrueInteraction&)> isSignal
  );

  const SpillMultiVar kNuMIPPFXWeightVectorPerSignalNu(
    [](const caf::SRSpillProxy *sr) -> vector<double> {
      return GetNuMIPPFXWeightVectorPerNu(sr, Is1muNp0piWithProtonPcut);
    }
  );

}
