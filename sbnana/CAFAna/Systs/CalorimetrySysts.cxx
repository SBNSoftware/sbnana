#include "sbnana/CAFAna/Systs/CalorimetrySysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <cmath>

namespace ana {

   CalorimetrySyst::CalorimetrySyst(CaloSyst _GainSyst, CaloSyst _AlphaSyst, CaloSyst _BetaSyst, const std::string& name, const std::string& latexName):
     ISyst(name, latexName),
     gain(0.01265), alpha(0.93), beta(0.212),
     gain_err(0.01), alpha_err(0.01), beta_err(0.01),
     GainSyst(_GainSyst), AlphaSyst(_AlphaSyst), BetaSyst(_BetaSyst)
  {
  }

  void CalorimetrySyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {

    double gain_new = (GainSyst==ana::CaloSyst::kNominal) ? gain : (1. + double(GainSyst)*gain_err)*gain;
    double alpha_new = (AlphaSyst==ana::CaloSyst::kNominal) ? alpha : (1. + double(AlphaSyst)*alpha_err)*alpha; // dummy
    double beta_new = (BetaSyst==ana::CaloSyst::kNominal) ? beta : (1. + double(BetaSyst)*beta_err)*beta;

    for(auto& pfp: sr->reco.pfp){
      for(int i_plane=0; i_plane<3; ++i_plane){
        // TODO ke
        for(auto& pt: pfp.trk.calo[i_plane].points){
          auto& this_dedx = pt.dedx;
          if(isnan(this_dedx)) continue;
          double new_dedx = std::exp( (beta_new*gain_new)/(beta*gain) * std::log(alpha + beta*this_dedx) - alpha_new  ) / beta_new;
          pt.dedx = new_dedx;
        }
      }
    }

  }

} // end namespace ana
