#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include <string>
#include <string_view>
#include <vector>

class TH1;

namespace ana
{

  class NuMIFluxSyst : public ISyst
  {
    static constexpr std::string_view fluxFileName = "2023-07-06_out_450.37_7991.98_79512.66_QEL11.root";

  public:
    virtual ~NuMIFluxSyst();

    void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;

  protected:
    friend const NuMIFluxSyst* GetNuMIFluxSyst(const std::string&,
                                               const std::string&,
                                               const std::string&);

    NuMIFluxSyst(const std::string& dir,
                             const std::string& prefix,
                             const std::string& name)
      : ISyst("numi_"+name, "NuMI flux: "+name),
        fHistName(dir+"/"+prefix+name), fName(name), fScale()
    {
    }

    std::string fHistName, fName;

    mutable TH1* fScale[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
  };

  const NuMIFluxSyst* GetNuMIFluxSyst(const std::string& dir,
                                      const std::string& prefix,
                                      const std::string& name);

  /// These are envelopes not real systs. TODO make clearer in naming
  std::vector<const ISyst*> GetNuMIHadronProductionFluxSysts();

  std::vector<const ISyst*> GetNuMIBeamlineFluxSysts();
  /// \param Npcs Number of principal components. High-indexed components
  ///             should have negligible effect on the analysis.
  std::vector<const ISyst*> GetNuMIPCAFluxSysts(unsigned int Npcs = 100);
  /// \brief Combination of all beamline systs plus \a Npcs hadron production
  ///        components
  std::vector<const ISyst*> GetAllNuMIFluxSysts(unsigned int Npcs);

} // namespace ana
