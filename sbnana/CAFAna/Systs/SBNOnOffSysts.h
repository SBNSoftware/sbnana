#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

#include <unordered_map>
#include <vector>

namespace ana
{
  class SBNOnOffSyst: public ISyst
  {
  public:
    SBNOnOffSyst(const std::string& systName);

    void Shift(double x, caf::SRSliceProxy* sr, double& weight) const override;

  protected:
    mutable int fIdx;
  };

  std::vector<std::string> GetSBNOnOffNames();
  const std::vector<const ISyst*>& GetSBNOnOffSysts();
}
