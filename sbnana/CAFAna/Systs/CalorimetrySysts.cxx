#include "sbnana/CAFAna/Systs/CalorimetrySysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "cetlib/search_path.h"
#include <cmath>
#include <iostream>
#include "TMath.h"

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

  Chi2Results CalorimetrySyst::CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo) const {

    // copied&modified from https://github.com/LArSoft/larana/blob/develop/larana/ParticleIdentification/Chi2PIDAlg.cxx#L60

    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double PIDA = 0; //by Bruce Baller
    std::vector<double> vpida;

    int used_trkres = 0;
    for(unsigned i = 0; i < calo.points.size(); ++i) { //hits
      const auto& pt = calo.points[i];
      double hit_dedx = pt.dedx;
      double hit_rr = pt.rr;

      //ignore the first and the last point
      if (i == 0 || i == calo.points.size() - 1) continue;
      if (hit_rr < 30) {
        PIDA += hit_dedx * pow(hit_rr, 0.42);
        vpida.push_back(hit_dedx * pow(hit_rr, 0.42));
        used_trkres++;
      }
      if (hit_dedx > 1000) continue; //protect against large pulse height
      int bin = dedx_range_pro->FindBin(hit_rr);
      if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
        double bincpro = dedx_range_pro->GetBinContent(bin);
        if (bincpro < 1e-6) { //for 0 bin content, using neighboring bins
          bincpro =
            (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
        }
        double bincka = dedx_range_ka->GetBinContent(bin);
        if (bincka < 1e-6) {
          bincka =
            (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;
        }
        double bincpi = dedx_range_pi->GetBinContent(bin);
        if (bincpi < 1e-6) {
          bincpi =
            (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
        }
        double bincmu = dedx_range_mu->GetBinContent(bin);
        if (bincmu < 1e-6) {
          bincmu =
            (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
        }
        double binepro = dedx_range_pro->GetBinError(bin);
        if (binepro < 1e-6) {
          binepro =
            (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;
        }
        double bineka = dedx_range_ka->GetBinError(bin);
        if (bineka < 1e-6) {
          bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;
        }
        double binepi = dedx_range_pi->GetBinError(bin);
        if (binepi < 1e-6) {
          binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;
        }
        double binemu = dedx_range_mu->GetBinError(bin);
        if (binemu < 1e-6) {
          binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;
        }
        //double errke = 0.05*hit_dedx;   //5% KE resolution
        double errdedx = 0.04231 + 0.0001783 * hit_dedx * hit_dedx; //resolution on dE/dx
        errdedx *= hit_dedx;
        chi2pro += pow((hit_dedx - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
        chi2ka += pow((hit_dedx - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
        chi2pi += pow((hit_dedx - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
        chi2mu += pow((hit_dedx - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
        //std::cout<<i<<" "<<hit_dedx<<" "<<hit_rr<<" "<<bincpro<<std::endl;
        ++npt;
      }
    }

    // Making output
    Chi2Results output;
    if (npt) {
      chi2pro /= npt;
      chi2ka /= npt;
      chi2pi /= npt;
      chi2mu /= npt;
    }
    static bool fUseMedian = true; // I think default fcl in ICARUS is true
    if (used_trkres > 0) {
      if (fUseMedian) {
        PIDA = TMath::Median(vpida.size(), &vpida[0]);
      }
      else { // use mean
        PIDA /= used_trkres;
      }
    }

    output.chi2_kaon = chi2ka;
    output.chi2_muon = chi2mu;
    output.chi2_pion = chi2pi;
    output.chi2_proton = chi2pro;
    output.pida = PIDA;
    output.pid_ndof = npt;

    return output;

  }

  void CalorimetrySyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {

    // TODO sigma is ambiguous for this kind of multi-param cases; not using it for now

    double gain_new = (1. + double(GainSyst)*gain_err)*gain;
    double alpha_new = (1. + double(AlphaSyst)*alpha_err)*alpha; // dummy
    double beta_new = (1. + double(BetaSyst)*beta_err)*beta;

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
        Chi2Results output = CalculateChi2(pfp.trk.calo[i_plane]);
        pfp.trk.chi2pid[i_plane].chi2_proton = output.chi2_proton;
        pfp.trk.chi2pid[i_plane].chi2_kaon = output.chi2_kaon;
        pfp.trk.chi2pid[i_plane].chi2_pion = output.chi2_pion;
        pfp.trk.chi2pid[i_plane].chi2_muon = output.chi2_muon;
        pfp.trk.chi2pid[i_plane].pida = output.pida;
        pfp.trk.chi2pid[i_plane].pid_ndof = output.pid_ndof;

/*
        // FOR DEBUGGING
        if( fabs(pfp.trk.chi2pid[i_plane].chi2_proton-276.916)<0.01 ){

          int ndof_old = pfp.trk.chi2pid[i_plane].pid_ndof;
          double chi2_proton_old = pfp.trk.chi2pid[i_plane].chi2_proton;
          double chi2_kaon_old = pfp.trk.chi2pid[i_plane].chi2_kaon;
          double chi2_pion_old = pfp.trk.chi2pid[i_plane].chi2_pion;
          double chi2_muon_old = pfp.trk.chi2pid[i_plane].chi2_muon;

          Chi2Results output = CalculateChi2(pfp.trk.calo[i_plane]);
          std::cout << "===============" << std::endl;
          std::cout << "sr->is_clear_cosmic = " << sr->is_clear_cosmic << std::endl;
          std::cout << "pfp ID = " << pfp.id << std::endl;
          std::cout << "plane = " << i_plane << std::endl;
          std::cout << "number of calo points = " << pfp.trk.calo[i_plane].points.size() << std::endl;
          printf("ndof, (old, new) = (%d, %d)\n",ndof_old,output.pid_ndof);
          printf("chi2_proton, (old, new) = (%1.3f, %1.3f)\n", chi2_proton_old, output.chi2_proton);
          printf("chi2_kaon, (old, new) = (%1.3f, %1.3f)\n", chi2_kaon_old, output.chi2_kaon);
          printf("chi2_pion, (old, new) = (%1.3f, %1.3f)\n", chi2_pion_old, output.chi2_pion);
          printf("chi2_muon, (old, new) = (%1.3f, %1.3f)\n", chi2_muon_old, output.chi2_muon);
        }
*/

      }
    }

  }

  void CalorimetrySyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

  }

  const CalorimetrySyst CalorimetrySyst_NoShift(ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, "CalorimetrySyst_NoShift", "NoShift");
  const CalorimetrySyst CalorimetrySyst_GainUp(ana::CaloSyst::kUp, ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, "CalorimetrySyst_GainUp", "gain +1%");
  const CalorimetrySyst CalorimetrySyst_GainDown(ana::CaloSyst::kDown, ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, "CalorimetrySyst_GainDown", "gain -1%");
  const CalorimetrySyst CalorimetrySyst_AlphaUp(ana::CaloSyst::kNominal, ana::CaloSyst::kUp, ana::CaloSyst::kNominal, "CalorimetrySyst_AlphaUp", "#alpha +4%");
  const CalorimetrySyst CalorimetrySyst_AlphaDown(ana::CaloSyst::kNominal, ana::CaloSyst::kDown, ana::CaloSyst::kNominal, "CalorimetrySyst_AlphaDown", "#alpha -4%");
  const CalorimetrySyst CalorimetrySyst_BetaUp(ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, ana::CaloSyst::kUp, "CalorimetrySyst_BetaUp", "#beta +4%");
  const CalorimetrySyst CalorimetrySyst_BetaDown(ana::CaloSyst::kNominal, ana::CaloSyst::kNominal, ana::CaloSyst::kDown, "CalorimetrySyst_BetaDown", "#beta -4%");
  const CalorimetrySyst CalorimetrySyst_GainUpBetaDown(ana::CaloSyst::kUp, ana::CaloSyst::kNominal, ana::CaloSyst::kDown, "CalorimetrySyst_GainUpBetaDown", "gain +1%, #beta -4%");
  const CalorimetrySyst CalorimetrySyst_GainDownBetaUp(ana::CaloSyst::kDown, ana::CaloSyst::kNominal, ana::CaloSyst::kUp, "CalorimetrySyst_GainDownBetaUp", "gain -1%, #beta +4%");

} // end namespace ana
