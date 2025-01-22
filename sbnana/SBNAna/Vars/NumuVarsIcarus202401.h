#pragma once

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TProfile.h"

#include <string>
#include <iostream>

namespace ana {

namespace chi2pid {

// Adapted from Jaesung's Chi2 Syst
struct Chi2PID {
  Chi2PID() : file_dEdXUncTemplate("template_dEdXUncertainty.root"), file_Chi2Template("dEdxrestemplates.root") {
    dedx_unc_template = (TGraph2D*)file_dEdXUncTemplate.Get("dEdXRelUncertainty_dEdX_vs_phi");
    dedx_range_pro = (TProfile*)file_Chi2Template.Get("dedx_range_pro");
    dedx_range_ka  = (TProfile*)file_Chi2Template.Get("dedx_range_ka");
    dedx_range_pi  = (TProfile*)file_Chi2Template.Get("dedx_range_pi");
    dedx_range_mu  = (TProfile*)file_Chi2Template.Get("dedx_range_mu");
  }

  struct Chi2PIDResult {
    double chi2_proton, chi2_kaon, chi2_pion, chi2_muon, PIDA;
    int ndof;
  };

  Chi2PIDResult calculate_chi2(const caf::Proxy<caf::SRTrackCalo> &calo) const {
    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double PIDA = 0; //by Bruce Baller
    std::vector<double> vpida;
    Chi2PIDResult output{0, 0, 0, 0, 0, 0};

    int used_trkres = 0;
    // If no calo points return all 0s
    if(calo.points.size() == 0) return output;
    //ignore the first and the last point
    for(unsigned i = 1; i < calo.points.size() - 1; ++i) { //hits
      const auto& pt = calo.points[i];
      double hit_dedx = pt.dedx;
      double hit_rr = pt.rr;

      // To match Maria, only consider last 25 cm
      if (hit_rr > 25) continue;
      if (hit_rr < 25) {
        PIDA += hit_dedx * pow(hit_rr, 0.42);
        vpida.push_back(hit_dedx * pow(hit_rr, 0.42));
        used_trkres++;
      }
      if (hit_dedx > 100 || hit_dedx < 0.5) continue; //protect against large pulse height
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

    // Making output
    output.chi2_kaon = chi2ka;
    output.chi2_muon = chi2mu;
    output.chi2_pion = chi2pi;
    output.chi2_proton = chi2pro;
    output.PIDA = PIDA;
    output.ndof = npt;

    return output;
  }

  TFile file_dEdXUncTemplate, file_Chi2Template;
  TGraph2D *dedx_unc_template;
  TProfile *dedx_range_pro, *dedx_range_ka, *dedx_range_pi, *dedx_range_mu;
};

extern Chi2PID chi2_calculator;

}

extern const Var kIcarus202401MuonIdx;
extern const Var kIcarus202401NumPions;
extern const Var kIcarus202401NumProtons;
extern const Var kIcarus202401NumShowers;

extern const Var kIcarus202401RecoMuonE;
extern const Var kIcarus202401RecoProtonKE;
extern const Var kIcarus202401RecoPionE;

extern const Var kIcarus202401RecoENu;

extern const Var kIcarus202401RecoMuonP;
extern const MultiVar kIcarus202401RecoProtonP;
extern const MultiVar kIcarus202401RecoPionP;


}


