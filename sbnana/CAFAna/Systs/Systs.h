#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana
{

class MECSyst: public ISyst
{
  public:
  MECSyst() : ISyst("mec", "100% MEC Systematic") {}
  void Shift(double sigma, caf::SRSliceProxy* sr, double& weight) const override
  {
    return this->Shift(sigma, &sr->truth, weight);
  }
  void Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const override
  {
    if(sigma < -1) sigma = -1;
    if(nu->genie_mode == 10) weight *= 1 + sigma;
  }
};

const MECSyst& GetMECSyst();

class POTSyst: public ISyst
{
  public:
  POTSyst() : ISyst("potsyst", "2% Systematic on POT") {}
  void Shift(double sigma, caf::SRSliceProxy* sr, double& weight) const override
  {
    return this->Shift(sigma, &sr->truth, weight);
  }
  void Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const override
  {
    weight *= (1.0 + 0.02 * sigma);
  }
};

const POTSyst& GetPOTSyst();

}
