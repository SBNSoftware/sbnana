#pragma once

#include <string>

namespace osc
{
  template<class T> class _IOscCalcAdjustable;
  typedef _IOscCalcAdjustable<double> IOscCalcAdjustable;
}

namespace ana
{
  /// Interface definition for fittable variables
  class IFitVar
  {
  public:
    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const = 0;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const = 0;
    virtual std::string ShortName() const = 0;
    virtual std::string LatexName() const = 0;
    virtual double Penalty(double /*val*/,
                           osc::IOscCalcAdjustable*) const {return 0;}
  };

  //----------------------------------------------------------------------

  /// Base class for variables with constraints. Apply penalty outside range
  class IConstrainedFitVar: public IFitVar
  {
  public:
    virtual double Penalty(double val, osc::IOscCalcAdjustable*) const;
    virtual double LowLimit() const = 0;
    virtual double HighLimit() const = 0;
  protected:
    double Clamp(double val) const;
  };
} // namespace
