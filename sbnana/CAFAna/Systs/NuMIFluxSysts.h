#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include <string>
#include <vector>

class TH1;

namespace ana {

  class NuMIFluxSyst : public ISyst {
  public:
    virtual ~NuMIFluxSyst();

    void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const override;

  protected:
    friend const NuMIFluxSyst* GetNuMIFluxSyst(const std::string&,
                                               const std::string&,
                                               const std::string&,
                                               const std::string&);

    NuMIFluxSyst(const std::string& dir,
                 const std::string& prefix,
                 const std::string& name,
                 const std::string& fluxFile = "");

    std::string fHistName, fName, fFluxFilePath;

    mutable TH1* fScale[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
  };

  const NuMIFluxSyst* GetNuMIFluxSyst(const std::string& dir,
                                      const std::string& prefix,
                                      const std::string& name,
                                      const std::string& fluxFile = "");

  /// These are envelopes not real systs. TODO make clearer in naming
  const ISyst* GetNuMIFluxStatisticalUncertainty();

  std::vector<const ISyst*> GetNuMIHadronProductionFluxSysts();

  std::vector<const ISyst*> GetNuMIBeamlineFluxSysts();
  /// \param Npcs Number of principal components. High-indexed components
  ///             should have negligible effect on the analysis.
  std::vector<const ISyst*> GetNuMIPCAFluxSysts(unsigned int Npcs = 100);
  /// \brief Combination of all beamline systs plus \a Npcs hadron production
  ///        components
  std::vector<const ISyst*> GetAllNuMIFluxSysts(unsigned int Npcs);

} // namespace ana
