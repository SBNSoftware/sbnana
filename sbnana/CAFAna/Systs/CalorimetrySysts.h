//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include <vector>

namespace ana
{

  enum class CaloSyst {
    kDown=-1,
    kNominal=0,
    kUp=+1,
  };

  class CalorimetrySyst: public ISyst
  {
  public:

    CalorimetrySyst(CaloSyst _GainSyst, CaloSyst _AlphaSyst, CaloSyst _BetaSyst, const std::string& name, const std::string& latexName);

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;

    // Update Nominal parameters
    inline void UpdateGain(double g, double gerr){ gain = g; gain_err = gerr; };
    inline void UpdateAlpha(double a, double aerr){ alpha = a; alpha_err = aerr; };
    inline void UpdateBeta(double b, double berr){ beta = b; beta_err = berr; };

  private:
    // nominal parameters
    double gain, alpha, beta;
    double gain_err, alpha_err, beta_err; // fractional uncertainty
    CaloSyst GainSyst, AlphaSyst, BetaSyst;
  };

}
