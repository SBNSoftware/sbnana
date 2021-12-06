#include "sbnana/CAFAna/Analysis/CutOptimizer.h"

#include "TH1.h"

#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  void MakeNMinusOneSpectra(SpectrumLoader& loader,
                            const Cut& sigcut,
                            const Cut& presel,
                            const std::vector<HistAxis>& axes,
                            const std::vector<double>& cut_pos,
                            std::vector<Spectrum>& sigs,
                            std::vector<Spectrum>& bkgs)
  {
    assert(cut_pos.size() == axes.size());

    sigs.reserve(axes.size());
    bkgs.reserve(axes.size());

    for(unsigned int i = 0; i < axes.size(); ++i){
      Cut nminusone = presel;
      for(unsigned j = 0; j < axes.size(); ++j){
        if(j == i) continue;
        // TODO support cuts with a < sign
        nminusone = nminusone && axes[j].GetVars()[0] > cut_pos[j];
      } // end for j

      sigs.emplace_back(loader, axes[i], nminusone && sigcut);
      bkgs.emplace_back(loader, axes[i], nminusone && !sigcut);
    } // end for i
  }

  // --------------------------------------------------------------------------
  double FindOptimumCut(TH1* hsig, TH1* hbkg, double& best_fom)
  {
    double nsig = 0;
    double nbkg = 0;
    best_fom = 0;
    double best_cut = -1;

    for(int binIdx = hsig->GetNbinsX()+1; binIdx >= 0; --binIdx){
      nsig += hsig->GetBinContent(binIdx);
      nbkg += hbkg->GetBinContent(binIdx);
      const double fom = nsig ? nsig/sqrt(nsig+nbkg) : 0;
      if(fom > best_fom){
        best_fom = fom;
        best_cut = hsig->GetXaxis()->GetBinLowEdge(binIdx);
      }
    } // end for binIdx

    return best_cut;
  }

  // --------------------------------------------------------------------------
  double OptimizeOneCut(const std::string& wildcard,
                        double pot,
                        const Cut& sigcut,
                        const Cut& presel,
                        const std::vector<HistAxis>& axes,
                        std::vector<double>& cut_pos)
  {
    SpectrumLoader loader(wildcard);
    std::vector<Spectrum> sigs, bkgs;
    MakeNMinusOneSpectra(loader, sigcut, presel, axes, cut_pos, sigs, bkgs);
    loader.Go();

    double best_fom = 0;
    double best_cut = -1;
    int best_idx = -1;

    for(unsigned int cutIdx = 0; cutIdx < axes.size(); ++cutIdx){
      TH1* hsig = sigs[cutIdx].ToTH1(pot, kRed);
      TH1* hbkg = bkgs[cutIdx].ToTH1(pot, kBlue);
      // TODO support multiple background spectra, such as cosmics

      double fom;
      const double cut = FindOptimumCut(hsig, hbkg, fom);

      if(fom > best_fom){
        best_fom = fom;
        best_cut = cut;
        best_idx = cutIdx;
      }

      // TODO plot distributions and optimized cuts at each phase
    } // end for cutIdx

    std::cout << "Updated cut on '" << axes[best_idx].GetLabels()[0] << "' to " << best_cut << ". FOM now " << best_fom << std::endl;

    // Store the result
    cut_pos[best_idx] = best_cut;
    return best_fom;
  }

  // --------------------------------------------------------------------------
  void OptimizeCuts(const std::string& wildcard,
                    double pot,
                    const Cut& sigcut,
                    const Cut& presel,
                    const std::vector<HistAxis>& axes,
                    std::vector<double>& cut_pos)
  {
    std::cout << "Initial cuts:" << std::endl;
    for(unsigned int i = 0; i < axes.size(); ++i){
      std::cout << "  " << axes[i].GetLabels()[0] << " > " << cut_pos[i] << std::endl;
    }

    double fom = 0;
    while(true){
      const double new_fom = OptimizeOneCut(wildcard, pot, sigcut, presel, axes, cut_pos);
      // Give up once the new FOM is not more than a 1% improvement
      if(new_fom < fom*1.01) break;
      fom = new_fom;
    }

    std::cout << "Final optimized cuts:" << std::endl;
    for(unsigned int i = 0; i < axes.size(); ++i){
      std::cout << "  " << axes[i].GetLabels()[0] << " > " << cut_pos[i] << std::endl;
    }
  }
}
