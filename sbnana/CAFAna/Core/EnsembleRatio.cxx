#include "sbnana/CAFAna/Core/EnsembleRatio.h"

#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

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
