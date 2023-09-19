#include "sbnana/CAFAna/Systs/CalorimetrySysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "cetlib/search_path.h"
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

    cet::search_path sp("FW_SEARCH_PATH");

    std::string fTemplateFile = "dEdxrestemplates.root";

    sp.find_file(fTemplateFile, fROOTfile);

    TFile *file = TFile::Open(fROOTfile.c_str());
    dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
    dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
    dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
    dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

  }

  double CalorimetrySyst::GetdEdXFromTemplate(double rr, ParticleType ptlType) const {

    int bin = dedx_range_pro->FindBin(rr);
    if(bin>dedx_range_pro->GetNbinsX()) bin = dedx_range_pro->GetNbinsX()-1;
    if(bin>=1&&bin<=dedx_range_pro->GetNbinsX()){

      double bincpro = dedx_range_pro->GetBinContent(bin);
      if (bincpro<1e-6){//for 0 bin content, using neighboring bins
        bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2;
      }
      double bincka = dedx_range_ka->GetBinContent(bin);
      if (bincka<1e-6){
        bincka = (dedx_range_ka->GetBinContent(bin-1)+dedx_range_ka->GetBinContent(bin+1))/2;
      }
      double bincpi = dedx_range_pi->GetBinContent(bin);
      if (bincpi<1e-6){
        bincpi = (dedx_range_pi->GetBinContent(bin-1)+dedx_range_pi->GetBinContent(bin+1))/2;
      }
      double bincmu = dedx_range_mu->GetBinContent(bin);
      if (bincmu<1e-6){
        bincmu = (dedx_range_mu->GetBinContent(bin-1)+dedx_range_mu->GetBinContent(bin+1))/2;
      }

      if(ptlType==ana::ParticleType::kProton) return bincpro;
      else if(ptlType==ana::ParticleType::kKaon) return bincka;
      else if(ptlType==ana::ParticleType::kPion) return bincpi;
      else if(ptlType==ana::ParticleType::kMuon) return bincmu;
      else{
        abort();
        return -1.;
      }

    }
    else{
      return -1.;
    }

  }

  double CalorimetrySyst::GetdEdXErrFromTemplate(double rr, ParticleType ptlType) const {

    int bin = dedx_range_pro->FindBin(rr);
    if(bin>dedx_range_pro->GetNbinsX()) bin = dedx_range_pro->GetNbinsX()-1;
    if(bin>=1&&bin<=dedx_range_pro->GetNbinsX()){

      double binepro = dedx_range_pro->GetBinError(bin);
      if (binepro<1e-6){
        binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2;
      }
      double bineka = dedx_range_ka->GetBinError(bin);
      if (bineka<1e-6){
        bineka = (dedx_range_ka->GetBinError(bin-1)+dedx_range_ka->GetBinError(bin+1))/2;
      }
      double binepi = dedx_range_pi->GetBinError(bin);
      if (binepi<1e-6){
        binepi = (dedx_range_pi->GetBinError(bin-1)+dedx_range_pi->GetBinError(bin+1))/2;
      }
      double binemu = dedx_range_mu->GetBinError(bin);
      if (binemu<1e-6){
        binemu = (dedx_range_mu->GetBinError(bin-1)+dedx_range_mu->GetBinError(bin+1))/2;
      }

      if(ptlType==ana::ParticleType::kProton) return binepro;
      else if(ptlType==ana::ParticleType::kKaon) return bineka;
      else if(ptlType==ana::ParticleType::kPion) return binepi;
      else if(ptlType==ana::ParticleType::kMuon) return binemu;
      else{
        abort();
        return -1.;
      }

    }
    else{
      return -1.;
    }

  }

  double CalorimetrySyst::CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo, ParticleType ptlType) const {

    int npt = 0;
    double chi2 = 0;
    for(unsigned i = 0; i < calo.points.size(); ++i) { //hits
      const auto& pt = calo.points[i];
      if (i == 0 || i == calo.points.size() - 1) continue;

      double dedx = GetdEdXFromTemplate(pt.rr, ptlType);
      double dedx_err = GetdEdXErrFromTemplate(pt.rr, ptlType);

      double errdedx = 0.04231 + 0.0001783 * pt.dedx * pt.dedx; //resolution on dE/dx
      errdedx *= pt.dedx;
      chi2 += pow( (pt.dedx-dedx)/std::sqrt( pow(dedx_err, 2) + pow(errdedx, 2) ), 2);
      npt++;

    }
    if(npt){
      chi2 /= npt;
    }

    return chi2;

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

        // shift dedx
        for(auto& pt: pfp.trk.calo[i_plane].points){
          auto& this_dedx = pt.dedx;
          if(isnan(this_dedx)) continue;
          if(isinf(this_dedx)) continue;
          double new_dedx = ( std::exp( (beta_new*gain_new)/(beta*gain) * std::log(alpha + (beta/rho/Efield)*this_dedx) ) - alpha_new ) / (beta_new/rho/Efield);
          pt.dedx = new_dedx;
        }
        // shift chi2
        pfp.trk.chi2pid[i_plane].chi2_proton = CalculateChi2(pfp.trk.calo[i_plane], ana::ParticleType::kProton);
        pfp.trk.chi2pid[i_plane].chi2_kaon = CalculateChi2(pfp.trk.calo[i_plane], ana::ParticleType::kKaon);
        pfp.trk.chi2pid[i_plane].chi2_pion = CalculateChi2(pfp.trk.calo[i_plane], ana::ParticleType::kPion);
        pfp.trk.chi2pid[i_plane].chi2_muon = CalculateChi2(pfp.trk.calo[i_plane], ana::ParticleType::kMuon);


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
