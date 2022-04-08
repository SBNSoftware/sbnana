#include "sbnana/CAFAna/Core/EnsembleRatio.h"

#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

#include "sbnana/CAFAna/Core/Multiverse.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "TGraphAsymmErrors.h"
#include "TH1.h"

namespace ana
{
  //----------------------------------------------------------------------
  EnsembleRatio::EnsembleRatio(const EnsembleSpectrum& num,
                               const EnsembleSpectrum& denom)
    : fMultiverse(&num.GetMultiverse()),
      fHist(num.fHist),
      fAxis(num.GetLabels(), num.GetBinnings())
  {
    CheckMultiverses(denom.GetMultiverse(), __func__);

    fHist.Divide(denom.fHist);
    // TODO TODO this show handle livetime too
    fHist.Scale(denom.POT()/num.POT());

    // This is clumsy, but the old histogram operation considered 0/0 = 0,
    // which is actually pretty useful (at least PredictionInterp relies on
    // this).
    for(int i = 0; i < fHist.GetNbinsX()+2; ++i){
      if(denom.fHist.GetBinContent(i) == 0){
        if(num.fHist.GetBinContent(i) == 0){
          fHist.SetBinContent(i, 0);
        }
        else{
          // Actual infinities break ROOT plotting
          fHist.SetBinContent(i, 1e100);
          // As in fact do mererly very large numbers
          // fHist.SetBinContent(i, std::numeric_limits<double>::max());
        }
      }
    }
  }

  //----------------------------------------------------------------------
  unsigned int EnsembleRatio::NUniverses() const
  {
    return fMultiverse->NUniv();
  }

  //----------------------------------------------------------------------
  Ratio EnsembleRatio::Universe(unsigned int univIdx) const
  {
    const int nbins = fAxis.GetBins1D().NBins()+2;

    if(fHist.HasStan()){
      // TODO
      abort();/*
      return Ratio(Eigen::ArrayXstan(fHist.GetEigenStan().segment(nbins*univIdx, nbins)),
                   fAxis.GetLabels(), fAxis.GetBinnings());
              */
    }
    else{
      return Ratio(Eigen::ArrayXd(fHist.GetEigen().segment(nbins*univIdx, nbins)),
                   fAxis.GetLabels(), fAxis.GetBinnings());
    }
  }

  //----------------------------------------------------------------------
  TGraphAsymmErrors* EnsembleRatio::ErrorBand() const
  {
    // TODO implementation without histograms

    std::unique_ptr<TH1D> hnom(Nominal().ToTH1());

    std::vector<std::unique_ptr<TH1D>> hunivs;
    hunivs.reserve(NUniverses()-1);
    for(unsigned int univIdx = 1; univIdx < NUniverses(); ++univIdx)
      hunivs.emplace_back(Universe(univIdx).ToTH1());

    TGraphAsymmErrors* g = new TGraphAsymmErrors;

    for(int binIdx = 0; binIdx < hnom->GetNbinsX()+2; ++binIdx){
      const double xnom = hnom->GetXaxis()->GetBinCenter(binIdx);
      const double ynom = hnom->GetBinContent(binIdx);
      g->SetPoint(binIdx, xnom, ynom);

      const double dx = hnom->GetXaxis()->GetBinWidth(binIdx);

      std::vector<double> ys;
      ys.reserve(hunivs.size());
      for(const std::unique_ptr<TH1D>& hu: hunivs){
        ys.push_back(hu->GetBinContent(binIdx));
      }

      // 1 sigma
      const double y0 = FindQuantile(.5-0.6827/2, ys);
      const double y1 = FindQuantile(.5+0.6827/2, ys);

      // It's theoretically possible for the central value to be outside the
      // error bands - clamp to zero in that case
      g->SetPointError(binIdx, dx/2, dx/2,
                       std::max(ynom-y0, 0.),
                       std::max(y1-ynom, 0.));
    } // end for binIdx

    return g;
  }

  //----------------------------------------------------------------------
  void EnsembleRatio::CheckMultiverses(const Multiverse& rhs,
                                       const std::string& func) const
  {
    if(&rhs == fMultiverse) return;

    std::cout << "EnsembleRatio::" << func << ": attempting to combine two spectra made with different multiverses: " << std::endl;
    std::cout << "  " << fMultiverse->ShortName() << std::endl;
    std::cout << "vs" << std::endl;
    std::cout << rhs.ShortName() << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  EnsembleRatio& EnsembleRatio::operator*=(const EnsembleRatio& rhs)
  {
    CheckMultiverses(rhs.GetMultiverse(), __func__);

    fHist.Multiply(rhs.fHist);

    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleRatio EnsembleRatio::operator*(const EnsembleRatio& rhs) const
  {
    EnsembleRatio ret = *this;
    ret *= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  EnsembleRatio& EnsembleRatio::operator/=(const EnsembleRatio& rhs)
  {
    CheckMultiverses(rhs.GetMultiverse(), __func__);

    fHist.Divide(rhs.fHist);

    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleRatio EnsembleRatio::operator/(const EnsembleRatio& rhs) const
  {
    EnsembleRatio ret = *this;
    ret /= rhs;
    return ret;
  }
}
