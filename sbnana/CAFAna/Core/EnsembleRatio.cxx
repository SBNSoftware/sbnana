#include "sbnana/CAFAna/Core/EnsembleRatio.h"

#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

#include "TGraphAsymmErrors.h"
#include "TH1.h"

namespace ana
{
  //----------------------------------------------------------------------
  EnsembleRatio::EnsembleRatio(const EnsembleSpectrum& num,
                               const EnsembleSpectrum& denom)
    : fNom(num.Nominal(), denom.Nominal())
  {
    assert(num.NUniverses() == denom.NUniverses());
    fUnivs.reserve(num.NUniverses());
    for(unsigned int i = 0; i < num.NUniverses(); ++i)
      fUnivs.emplace_back(num.Universe(i), denom.Universe(i));
  }

  //----------------------------------------------------------------------
  TGraphAsymmErrors* EnsembleRatio::ErrorBand() const
  {
    std::unique_ptr<TH1D> hnom(fNom.ToTH1());

    std::vector<std::unique_ptr<TH1D>> hunivs;
    hunivs.reserve(fUnivs.size());
    for(const Ratio& u: fUnivs) hunivs.emplace_back(u.ToTH1());

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
  EnsembleRatio& EnsembleRatio::operator*=(const EnsembleRatio& rhs)
  {
    assert(rhs.NUniverses() == fUnivs.size());
    fNom *= rhs.fNom;
    for(unsigned int i = 0; i < fUnivs.size(); ++i)
      fUnivs[i] *= rhs.fUnivs[i];
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
    assert(rhs.NUniverses() == fUnivs.size());
    fNom /= rhs.fNom;
    for(unsigned int i = 0; i < fUnivs.size(); ++i)
      fUnivs[i] /= rhs.fUnivs[i];
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
