#pragma once

#include "CAFAna/Core/INamed.h"

namespace osc
{
  template<class T> class _IOscCalcAdjustable;
  typedef _IOscCalcAdjustable<double> IOscCalcAdjustable;
}

namespace ana
{
  /// Interface definition for fittable variables
  class IFitVar: public INamed
  {
  public:
    IFitVar(const std::string& shortName, const std::string& latexName);

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const = 0;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const = 0;
    virtual double Penalty(double /*val*/,
                           osc::IOscCalcAdjustable*) const {return 0;}
  };

  //----------------------------------------------------------------------

  /// Base class for variables with constraints. Apply penalty outside range
  class IConstrainedFitVar: public IFitVar
  {
  public:
    IConstrainedFitVar(const std::string& shortName, const std::string& latexName);
    virtual double Penalty(double val, osc::IOscCalcAdjustable*) const;
    virtual double LowLimit() const = 0;
    virtual double HighLimit() const = 0;
  protected:
    double Clamp(double val) const;
  };
} // namespace
