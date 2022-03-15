#pragma once

#include "sbnana/CAFAna/Core/Ratio.h"

#include "TGraphAsymmErrors.h"

#include <cassert>

namespace ana
{
  class EnsembleSpectrum;

  class EnsembleRatio
  {
  public:
    EnsembleRatio(const EnsembleSpectrum& nom, const EnsembleSpectrum& denom);

    Ratio Nominal() const {return fNom;}
    unsigned int NUniverses() const {return fUnivs.size();}
    Ratio Universe(unsigned int i) const
    {
      assert(i < fUnivs.size());
      return fUnivs[i];
    }

    TGraphAsymmErrors* ErrorBand() const;

    EnsembleRatio& operator*=(const EnsembleRatio& rhs);
    EnsembleRatio operator*(const EnsembleRatio& rhs) const;

    EnsembleRatio& operator/=(const EnsembleRatio& rhs);
    EnsembleRatio operator/(const EnsembleRatio& rhs) const;

  protected:
    Ratio fNom;
    std::vector<Ratio> fUnivs;
  };

  inline EnsembleRatio operator/(const EnsembleSpectrum& lhs,
                                 const EnsembleSpectrum& rhs)
  {
    return EnsembleRatio(lhs, rhs);
  }
}
