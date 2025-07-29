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

  class NuMIBeamG3ChaseSyst: public ISyst
  {
  public:

    NuMIBeamG3ChaseSyst(const std::string& name, const std::string& latexName);

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:

  };

  /// Class *like* NuMIFluxSyst but for the BeamShift systematic
  // 05/19/25 JK) Updated to use 2025-04-08_out_450.37_7991.98_79512.66.root;
  //              The CV in this file has beam width of 1.5 mm, while 2023NuMI Reprocessing
  //              file had CV 1.4mm. Here we want to get beam focusing uncertainties around 1.5mm
  //              Users of 2023NuMI reprocessing need to apply a CV correction of 1.5mm/1.4mm
  class NuMIBeamShiftSyst : public ISyst {
  public:
    virtual ~NuMIBeamShiftSyst();

    void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const override;

  protected:
    friend const NuMIBeamShiftSyst* GetNuMIBeamShiftSyst(const std::string&, const std::string&);

    NuMIBeamShiftSyst(const std::string& name,
                      const std::string& fluxFile = "");

    // For the beam focusing systematics
    // "fName" is the dial name
    std::string fName, fFluxFilePath;

    mutable TH1* fScale[2][2][2][2]; // [+sigma/-sigma][fhc/rhc][nue/numu][nu/nubar]
  };

  const NuMIBeamShiftSyst* GetNuMIBeamShiftSyst(const std::string& name,
                                                const std::string& fluxFile);
  std::vector<const ISyst*> GetNuMIBeamShiftSysts();

} // namespace ana
