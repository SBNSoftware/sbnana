#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Weight.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

#include <unordered_map>
#include <vector>

namespace ana
{
  class UniverseWeight
  {
  public:
    UniverseWeight(const std::vector<std::string>& systs, int univIdx);
    UniverseWeight(const std::vector<const ISyst*>& systs, int univIdx);

    double operator()(const caf::SRSliceProxy* sr) const;

  protected:
    std::vector<std::string> fNames;
    int fUnivIdx;
    mutable std::vector<unsigned int> fSystIdxs;
    mutable std::vector<int> fUnivOffsets;
    mutable std::vector<Cut> fUnivCuts;
  };

  Weight GetUniverseWeight(const std::string& syst, int univIdx)
  {
    return Weight(UniverseWeight(std::vector<std::string>(1, syst), univIdx));
  }

  Weight GetUniverseWeight(const std::vector<std::string>& systs, int univIdx)
  {
    return Weight(UniverseWeight(systs, univIdx));
  }

  Weight GetUniverseWeight(const std::vector<const ISyst*> systs, int univIdx)
  {
    return Weight(UniverseWeight(systs, univIdx));
  }


  class SBNWeightSyst: public ISyst
  {
  public:
    SBNWeightSyst(const std::string& systName,
                  const std::string& knobName = "", // if it differs
                  const SliceCut& cut = kNoCut);

    void Shift(double x, caf::SRSliceProxy* sr, double& weight) const override;

  protected:
    mutable int fIdx;

    std::string fKnobName;

    SliceCut fCut;

    struct Univs
    {
      int i0, i1;
      double w0, w1;
    };

    mutable std::unordered_map<double, Univs> fUnivs;

    Univs GetUnivs(double x) const;
  };

  const std::vector<const ISyst*>& GetSBNGenieWeightSysts();

  //  const std::vector<const ISyst*>& GetSBNBoosterWeightSysts();
  //  const std::vector<const ISyst*>& GetSBNWeightSysts(); // genie+flux
}
