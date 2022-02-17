#pragma once

#include "sbnana/CAFAna/Core/IFitVar.h"

#pragma once

namespace ana
{
  /// \f$ \Delta m^2 \f$
  class FitDmSqSterile: public IConstrainedFitVar
  {
  public:
    FitDmSqSterile() : IConstrainedFitVar("dmsq", "#Deltam^{2} (eV^{2})") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1e6;}
  };

  /// \Delta m^2 \f$
  const FitDmSqSterile kFitDmSqSterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{\mu\mu} \f$
  class FitSinSq2ThetaMuMu: public IConstrainedFitVar
  {
  public:
    FitSinSq2ThetaMuMu() : IConstrainedFitVar("ss2thmm", "sin^{2}2#theta_{#mu#mu}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{\mu\mu} \f$
  const FitSinSq2ThetaMuMu kFitSinSq2ThetaMuMu;

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{\mu e} \f$
  class FitSinSq2ThetaMuE: public IConstrainedFitVar
  {
  public:
    FitSinSq2ThetaMuE() : IConstrainedFitVar("ss2thme", "sin^{2}2#theta_{#mue}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{\mu e} \f$
  const FitSinSq2ThetaMuE kFitSinSq2ThetaMuE;

} // namespace
