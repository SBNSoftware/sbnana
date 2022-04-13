#pragma once

#include "cafanacore/Ratio.h"

#include "TGraphAsymmErrors.h"

#include <cassert>

namespace ana
{
  class EnsembleSpectrum;
  class FitMultiverse;

  class EnsembleRatio
  {
  public:
    friend class EnsembleSpectrum;

    EnsembleRatio(const EnsembleSpectrum& num, const EnsembleSpectrum& denom);

    Ratio Nominal() const {return Universe(0);}
    unsigned int NUniverses() const;
    Ratio Universe(unsigned int i) const;

    const FitMultiverse& GetMultiverse() const {return *fMultiverse;}

    TGraphAsymmErrors* ErrorBand() const;

    EnsembleRatio& operator*=(const EnsembleRatio& rhs);
    EnsembleRatio operator*(const EnsembleRatio& rhs) const;

    EnsembleRatio& operator/=(const EnsembleRatio& rhs);
    EnsembleRatio operator/(const EnsembleRatio& rhs) const;

  protected:
    void CheckMultiverses(const FitMultiverse& rhs,
                          const std::string& func) const;

    const FitMultiverse* fMultiverse;
    Hist fHist;
    LabelsAndBins fAxis;
  };

  inline EnsembleRatio operator/(const EnsembleSpectrum& lhs,
                                 const EnsembleSpectrum& rhs)
  {
    return EnsembleRatio(lhs, rhs);
  }
}
