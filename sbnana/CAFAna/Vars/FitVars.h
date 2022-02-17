#pragma once

#include "sbnana/CAFAna/Core/IFitVar.h"

#include <limits>

namespace ana
{
  /// \f$ \theta_{13} \f$
  class FitTheta13: public IFitVar
  {
  public:
    FitTheta13() : IFitVar("th13", "#theta_{13}") {}
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \theta_{13} \f$
  const FitTheta13 kFitTheta13;

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{13} \f$
  class FitSinSq2Theta13: public IConstrainedFitVar
  {
  public:
    FitSinSq2Theta13() : IConstrainedFitVar("ss2th13", "sin^{2}2#theta_{13}") {}
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{13} \f$
  const FitSinSq2Theta13 kFitSinSq2Theta13;

  //----------------------------------------------------------------------

  /// \f$ \delta_{CP}/\pi \f$
  class FitDeltaInPiUnits: public IFitVar
  {
  public:
    FitDeltaInPiUnits() : IFitVar("delta(pi)", "#delta / #pi") {}
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \delta_{CP}/\pi \f$
  const FitDeltaInPiUnits kFitDeltaInPiUnits;

  //----------------------------------------------------------------------
  /// \f$ \theta_{13} \f$
  class FitTheta23: public IFitVar
  {
  public:
    FitTheta23() : IFitVar("th23", "#theta_{23}") {}
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \theta_{13} \f$
  const FitTheta23 kFitTheta23;
  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{23} \f$
  class FitSinSqTheta23: public IConstrainedFitVar
  {
  public:
    FitSinSqTheta23() : IConstrainedFitVar("ssth23", "sin^{2}#theta_{23}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{23} \f$
  const FitSinSqTheta23 kFitSinSqTheta23;

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{23} \f$
  class FitSinSq2Theta23: public IConstrainedFitVar
  {
  public:
    FitSinSq2Theta23() : IConstrainedFitVar("ss2th23", "sin^{2}2#theta_{23}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{23} \f$
  const FitSinSq2Theta23 kFitSinSq2Theta23;

  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{32} \f$
  class FitDmSq32: public IConstrainedFitVar
  {
  public:
    FitDmSq32() : IConstrainedFitVar("dmsq32", "#Deltam^{2}_{32}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    // "1eV^2 splitting should be enough for anyone"
    // OscCalcPMNS freaks out at large splittings
    virtual double LowLimit() const {return -1;}
    virtual double HighLimit() const {return +1;}
  };

  /// \f$ \Delta m^2_{32} \f$
  const FitDmSq32 kFitDmSq32;

  //-------------------------------------------------------------------------

  /// \f$ \Delta m^2_{32}\times10^3{\rm eV}^2 \f$
  class FitDmSq32Scaled: public IConstrainedFitVar
  {
  public:
    FitDmSq32Scaled() : IConstrainedFitVar("dmsq32scaled", "#Deltam^{2}_{32} (10^{-3} eV^{2})") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    // "1eV^2 splitting should be enough for anyone"
    // OscCalcPMNS freaks out at large splittings
    virtual double LowLimit() const {return -1000;}
    virtual double HighLimit() const {return +1000;}
  };

  /// \f$ \Delta m^2_{32}\times10^3{\rm eV}^2 \f$
  const FitDmSq32Scaled kFitDmSq32Scaled;

  //----------------------------------------------------------------------

  /// \f$ \tan^2\theta_{12} \f$
  class FitTanSqTheta12: public IConstrainedFitVar
  {
  public:
    FitTanSqTheta12() : IConstrainedFitVar("tsth12", "tan^{2}#theta_{12}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return std::numeric_limits<double>::max();}
  };

  /// \f$ \tan^2\theta_{12} \f$
  const FitTanSqTheta12 kFitTanSqTheta12;

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{12} \f$
  class FitSinSq2Theta12: public IConstrainedFitVar
  {
  public:
    FitSinSq2Theta12() : IConstrainedFitVar("ss2th12", "sin^{2}2#theta_{12}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{12} \f$
  const FitSinSq2Theta12 kFitSinSq2Theta12;

  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{21} \f$
  class FitDmSq21: public IConstrainedFitVar
  {
  public:
    FitDmSq21() : IConstrainedFitVar("dmsq21","#Deltam^{2}_{21}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    // "1eV^2 splitting should be enough for anyone"
    // OscCalcPMNS freaks out at large splittings
    virtual double LowLimit() const {return -1;}
    virtual double HighLimit() const {return +1;}
  };

  /// \f$ \Delta m^2_{21} \f$
  const FitDmSq21 kFitDmSq21;

  /// \f$ \rho \f$
  class FitRho: public IConstrainedFitVar
  {
  public:
    FitRho() : IConstrainedFitVar("rho", "#rho") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    //Density should be greater than zero (set a ridiculously high high limit)
    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 10.0;}

  };

  /// \f$ \rho \f$
  const FitRho kFitRho;

  //----------------------------------------------------------------------



} // namespace
