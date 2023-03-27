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
    UniverseWeight(const std::string& psetName, double x);

    double operator()(const caf::SRSpillProxy* sr) const;
    double operator()(const caf::SRSliceProxy* sr) const;

  protected:
    std::string fPSetName;
    mutable int fPSetIdx;
    int fUnivIdx;
    double fSigma;
  };

  SpillVar GetUniverseFirstNeutrinoWeight(const std::string& psetName, double x)
  {
    return SpillVar(UniverseWeight(psetName, x));
  }

  SpillVar GetUniverseWeightSpill(const std::string& psetName, int univIdx)
  {
    return SpillVar(UniverseWeight(psetName, univIdx));
  }
  Var GetUniverseWeight(const std::string& psetName, int univIdx)
  {
    return Var(UniverseWeight(psetName, univIdx));
  }


  class SBNWeightSyst: public ISyst
  {
  public:
    SBNWeightSyst(const std::string& systName);

    void Shift(double x, caf::SRSliceProxy* sr, double& weight) const override;

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
  std::vector<std::string> GetSBNBoosterFluxWeightNames();

  const std::vector<const ISyst*>& GetSBNGenieWeightSysts();

  const std::vector<const ISyst*>& GetSBNBoosterFluxWeightSysts();
  const std::vector<const ISyst*>& GetSBNWeightSysts(); // genie+flux
}
