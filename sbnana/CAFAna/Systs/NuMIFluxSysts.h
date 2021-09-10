#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include <vector>

class TH1;

namespace ana
{

  class NuMIFluxSyst : public ISyst
  {
  public:
    virtual ~NuMIFluxSyst();

    void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;

  protected:
    friend const NuMIFluxSyst* GetNuMIFluxSyst(const std::string&,
                                               const std::string&);

    NuMIFluxSyst(const std::string& dir, const std::string& name)
      : ISyst("numi_"+name, "NuMI flux: "+name),
        fDir(dir), fName(name), fScale()
    {
    }

    std::string fDir, fName;

    mutable TH1* fScale[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
  };

  const NuMIFluxSyst* GetNuMIFluxSyst(const std::string& dir,
                                      const std::string& name);

  std::vector<const ISyst*> GetNuMIHadronProductionFluxSysts();
  std::vector<const ISyst*> GetNuMIBeamlineFluxSysts();
  std::vector<const ISyst*> GetAllNuMIFluxSysts();

} // namespace ana
