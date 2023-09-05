#pragma once

#include "sbnana/CAFAna/Systs/UniverseOracle.h"

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

namespace ana{

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
