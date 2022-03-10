#pragma once

#include "cafanacore/ISyst.h"

#include <string>

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  /// \brief Encapsulate code to systematically shift a \ref
  /// caf::StandardRecord
  ///
  /// The Shift() function alters the \ref caf::StandardRecord or the weight
  /// associated with the event.
  class ISyst: public _ISyst<caf::SRSliceProxy>
  {
  public:
    ISyst(const std::string& shortName,
          const std::string& latexName,
	  bool applyPenalty = true,
	  double min = -3,
	  double max = +3);
    virtual ~ISyst();

    virtual double Penalty(double x) const;

    /// Should a penalty be applied for this shift?
    virtual bool ApplyPenalty() const {return fApplyPenalty;}

    /// Return the min/max value for this syst
    virtual double Min() const{return fMin;}
    virtual double Max() const{return fMax;}

    /// PredictionInterp normally interpolates between spectra made at
    /// +/-1,2,3sigma. For some systematics that's overkill. Override this
    /// function to specify different behaviour for this systematic.
    virtual int PredInterpMaxNSigma() const
    {
      return 3;
    }

  private:
    bool fApplyPenalty;
    double fMin;
    double fMax;
  };

} // namespace
