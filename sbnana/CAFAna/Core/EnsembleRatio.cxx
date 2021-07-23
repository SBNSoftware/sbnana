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

    TGraphAsymmErrors* g = new TGraphAsymmErrors;

    for(int binIdx = 0; binIdx < hnom->GetNbinsX()+2; ++binIdx){
      const double xnom = hnom->GetXaxis()->GetBinCenter(binIdx);
      const double ynom = hnom->GetBinContent(binIdx);
      g->SetPoint(binIdx, xnom, ynom);

      const double dx = hnom->GetXaxis()->GetBinWidth(binIdx);

      g->SetPointError(binIdx, dx/2, dx/2, .1*ynom, .1*ynom);
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
