//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TProfile.h"

#include <optional>

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
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:

    // detector parameters
    double temperature, rho, Efield;

    // nominal parameters
    double gain, alpha, beta;
    double gain_err;
    CaloSystMode kCaloSystMode;

    /*
     * Loading of the information is lazy and happens the first time
     * `rangeTemplates()` and `dEdXuncTemplates()` are called.
     * 
     */
    
    // dEdX uncertainty template
    struct dEdX_uncertainty_t {
      TGraph2D *templ = nullptr;
    };
    
    struct dEdX_range_templates_t {
      TProfile *pro = nullptr;    ///< proton template
      TProfile *ka  = nullptr;    ///< kaon template
      TProfile *pi  = nullptr;    ///< pion template
      TProfile *mu  = nullptr;    ///< muon template
    };

    /// Fetches and stores template information.
    class cache_t {
      
      mutable std::optional<dEdX_uncertainty_t> dedx_unc;       ///< dE/dX uncertainty
      mutable std::optional<dEdX_range_templates_t> dedx_range; ///< PID uncertainty templates
      
      /// Loads and returns the dE/dX uncertainty template(s).
      dEdX_uncertainty_t load_dEdXuncTemplates() const;
      
      /// Loads and returns the range templates.
      dEdX_range_templates_t load_rangeTemplates() const;
      
        public:
          
      /// Returns the range templates, loading them if needed.
      dEdX_uncertainty_t const& dEdXuncTemplates() const;
      
      /// Returns the range templates, loading them if needed.
      dEdX_range_templates_t const& rangeTemplates() const;
      
    }; // cache_t
    
    cache_t cache; ///< Data cache.
    
  };

  extern const CalorimetrySyst kCalodEdXShiftSyst;
  extern const CalorimetrySyst kCaloGainShiftSyst;

}
