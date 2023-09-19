#include "sbnana/CAFAna/Systs/CalorimetrySysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <cmath>

namespace ana {

   CalorimetrySyst::CalorimetrySyst(CaloSyst _GainSyst, CaloSyst _AlphaSyst, CaloSyst _BetaSyst, const std::string& name, const std::string& latexName):
     ISyst(name, latexName),
     temperature(8.75e1),
     rho(-0.00615 * temperature + 1.928),
     Efield(4.938e-1),
     gain(0.01265), alpha(0.93), beta(0.212),
     gain_err(0.01), alpha_err(0.01), beta_err(0.04),
     GainSyst(_GainSyst), AlphaSyst(_AlphaSyst), BetaSyst(_BetaSyst)
  {
  }

  void CalorimetrySyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {

    // TODO sigma is ambiguous for this kind of multi-param cases; not using it for now

    double gain_new = (GainSyst==ana::CaloSyst::kNominal) ? gain : (1. + double(GainSyst)*gain_err)*gain;
    double alpha_new = (AlphaSyst==ana::CaloSyst::kNominal) ? alpha : (1. + double(AlphaSyst)*alpha_err)*alpha; // dummy
    double beta_new = (BetaSyst==ana::CaloSyst::kNominal) ? beta : (1. + double(BetaSyst)*beta_err)*beta;

    for(auto& pfp: sr->reco.pfp){
      for(int i_plane=0; i_plane<3; ++i_plane){
        // TODO ke
        for(auto& pt: pfp.trk.calo[i_plane].points){
          auto& this_dedx = pt.dedx;
          if(isnan(this_dedx)) continue;
          if(isinf(this_dedx)) continue;
          double new_dedx = ( std::exp( (beta_new*gain_new)/(beta*gain) * std::log(alpha + (beta/rho/Efield)*this_dedx) ) - alpha_new ) / (beta_new/rho/Efield);
          pt.dedx = new_dedx;
        }
      }
    }

  }

  void CalorimetrySyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

  }

  const CalorimetrySyst CalorimetrySyst_BetaUp(ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, ana::CaloSyst::kUp, "CalorimetrySyst_BetaUp", "#beta +4%");
  const CalorimetrySyst CalorimetrySyst_BetaDown(ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, ana::CaloSyst::kDown, "CalorimetrySyst_BetaDown", "#beta -4%");

  const CalorimetrySyst CalorimetrySyst_GainUp(ana::CaloSyst::kUp, ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, "CalorimetrySyst_GainUp", "gain +1%");
  const CalorimetrySyst CalorimetrySyst_GainDown(ana::CaloSyst::kDown, ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, "CalorimetrySyst_GainUp", "gain +1%");

} // end namespace ana
