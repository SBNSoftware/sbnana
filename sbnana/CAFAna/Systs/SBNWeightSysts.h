#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

#include <unordered_map>
#include <vector>

namespace ana
{
  class UniverseWeight
  {
  public:
    UniverseWeight(const std::string& psetName, int univIdx);

    double operator()(const caf::SRSliceProxy* sr) const;
    double operator()(const caf::SRTrueInteractionProxy* nu) const;

  protected:
    std::string fPSetName;
    mutable int fPSetIdx;
    int fUnivIdx;
  };

  Var GetUniverseWeight(const std::string& psetName, int univIdx)
  {
    return Var(UniverseWeight(psetName, univIdx));
  }

  TruthVar GetTruthUniverseWeight(const std::string& psetName, int univIdx)
  {
    return TruthVar(UniverseWeight(psetName, univIdx));
  }


  class SBNWeightSyst: public ISyst
  {
  public:
    SBNWeightSyst(const std::string& systName);

    void Shift(double x, caf::SRSliceProxy* sr, double& weight) const override;
    void Shift(double x, caf::SRTrueInteractionProxy* sr, double& weight) const override;

  protected:
    mutable int fIdx;

    struct Univs
    {
      int i0, i1;
      double w0, w1;
    };

    mutable std::unordered_map<double, Univs> fUnivs;

    Univs GetUnivs(double x) const;
  };

  std::vector<std::string> GetSBNGenieWeightNames();

  const std::vector<const ISyst*>& GetSBNGenieWeightSysts();

  const std::vector<const ISyst*>& GetSBNBoosterFluxWeightSysts();
  const std::vector<const ISyst*>& GetSBNWeightSysts(); // genie+flux
}
