//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "TProfile.h"
#include "TFile.h"

#include <vector>

namespace ana
{

  enum CaloSyst {
    kDown=-1,
    kNominal=0,
    kUp=+1,
  };

  struct Chi2Results { ///< determined particle ID
    Float_t chi2_kaon, chi2_muon, chi2_pion, chi2_proton, pida;
    Int_t pid_ndof;
  };

  class CalorimetrySyst: public ISyst
  {
  public:

    CalorimetrySyst(CaloSyst _GainSyst, CaloSyst _AlphaSyst, CaloSyst _BetaSyst, const std::string& name, const std::string& latexName);

    Chi2Results CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo) const;

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

    // Update detector parameters
    inline void UpdateTemperature(double t){ temperature = t; }
    inline void UpdateDensity(double r){ rho = r; }
    inline void UpdateEfield(double e){ Efield = e; }

    // Update nominal parameters
    inline void UpdateGain(double g, double gerr){ gain = g; gain_err = gerr; }
    inline void UpdateAlpha(double a, double aerr){ alpha = a; alpha_err = aerr; }
    inline void UpdateBeta(double b, double berr){ beta = b; beta_err = berr; }

  private:

    // detector parameters
    double temperature, rho, Efield;

    // nominal parameters
    double gain, alpha, beta;
    double gain_err, alpha_err, beta_err; // fractional uncertainty
    CaloSyst GainSyst, AlphaSyst, BetaSyst;

    // chi2 calculation
    std::string fTemplateFile;
    std::string fROOTfile;

    TProfile *dedx_range_pro;   ///< proton template
    TProfile *dedx_range_ka;    ///< kaon template
    TProfile *dedx_range_pi;    ///< pion template
    TProfile *dedx_range_mu;    ///< muon template

  };

  extern const CalorimetrySyst CalorimetrySyst_NoShift; // for debugging
  extern const CalorimetrySyst CalorimetrySyst_GainUp;
  extern const CalorimetrySyst CalorimetrySyst_GainDown;
  extern const CalorimetrySyst CalorimetrySyst_AlphaUp;
  extern const CalorimetrySyst CalorimetrySyst_AlphaDown;
  extern const CalorimetrySyst CalorimetrySyst_BetaUp;
  extern const CalorimetrySyst CalorimetrySyst_BetaDown;
  extern const CalorimetrySyst CalorimetrySyst_GainUpBetaDown;
  extern const CalorimetrySyst CalorimetrySyst_GainDownBetaUp;

}
