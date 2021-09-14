#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include <vector>

namespace ana
{
  class MECSyst: public ISyst
  {
  public:
    MECSyst() : ISyst("mec", "MEC Syst") {}
    void Shift(double sigma, caf::SRSliceProxy* sr, double& weight) const override
    {
      if(sigma < -1) sigma = -1.0;
      if(sr->truth.genie_intcode == 10) weight *= 1.0 + sigma;
    }
  };

  const MECSyst& GetMECSyst()
  {
    static const MECSyst mSyst;
    return mSyst;
  }

  class POTSyst: public ISyst
  {
  public:
    POTSyst() : ISyst("potsyst", "2% POT Syst") {}
    void Shift(double sigma, caf::SRSliceProxy* sr, double& weight) const override
    {
      weight *= (1.0 + 0.02*sigma);
    }
  };

  const POTSyst& GetPOTSyst()
  {
    static const POTSyst pSyst;
    return pSyst;
  }

  class NormalizationSyst: public ISyst
  {
  public:
    NormalizationSyst() : ISyst("normsyst", "1% Normalization Syst") {}
    void Shift(double sigma, caf::SRSliceProxy *sr, double &weight) const override
    {
      weight *= (1.0 + 0.01*sigma);
    }
  };

  const NormalizationSyst& GetNormSyst()
  {
    static const NormalizationSyst nSyst;
    return nSyst;
  }

  class NormalizationSystND: public ISyst
  {
  public:
    NormalizationSystND() : ISyst("normsystnd", "1% Normalization Syst in SBND") {}
    void Shift(double sigma, caf::SRSliceProxy *sr, double &weight) const override
    {
      if(!isnan(sr->truth.baseline) && sr->truth.baseline < 120)
        weight *= (1.0 + 0.01*sigma);
    }
  };

  const NormalizationSystND& GetNormSystND()
  {
    static const NormalizationSystND nSyst;
    return nSyst;
  }

  class NormalizationSystFD: public ISyst
  {
  public:
    NormalizationSystFD() : ISyst("normsystfd", "1% Normalization Syst in ICARUS") {}
    void Shift(double sigma, caf::SRSliceProxy *sr, double &weight) const override
    {
      if(!isnan(sr->truth.baseline) && sr->truth.baseline > 500)
        weight *= (1.0 + 0.01*sigma);
    }
  };

  const NormalizationSystFD& GetNormSystFD()
  {
    static const NormalizationSystFD nSyst;
    return nSyst;
  }

  const std::vector<const ISyst*>& GetMiscSysts();
}
