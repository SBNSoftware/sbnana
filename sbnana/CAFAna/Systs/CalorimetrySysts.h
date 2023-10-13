//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TProfile.h"
#include "TMath.h"

#include <vector>

namespace ana
{

  enum CaloSystMode {
    kdEdXShift=0,
    kGainShift=1,
  };

  struct Chi2Results { ///< determined particle ID
    Float_t chi2_kaon, chi2_muon, chi2_pion, chi2_proton, pida;
    Int_t pid_ndof;
  };

  class CalorimetrySyst: public ISyst
  {
  public:

    CalorimetrySyst(CaloSystMode mode, const std::string& name, const std::string& latexName);

    Chi2Results CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo) const;

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;

  private:

    // detector parameters
    double temperature, rho, Efield;

    // nominal parameters
    double gain, alpha, beta;
    double gain_err;
    CaloSystMode kCaloSystMode;

    // dEdX uncertainty template
    TGraph2D *dedx_unc_template;

    TProfile *dedx_range_pro;   ///< proton template
    TProfile *dedx_range_ka;    ///< kaon template
    TProfile *dedx_range_pi;    ///< pion template
    TProfile *dedx_range_mu;    ///< muon template

  };

  extern const CalorimetrySyst kCalodEdXShiftSyst;
  extern const CalorimetrySyst kCaloGainShiftSyst;

}
