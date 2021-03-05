#include "sbnana/CAFAna/Experiment/GaussianConstraint.h"

#include "sbnana/CAFAna/Vars/FitVars.h"

#include "sbnana/CAFAna/Core/MathUtil.h"

namespace ana
{
  //----------------------------------------------------------------------
  double GaussianConstraint::ChiSq(osc::IOscCalcAdjustable* osc,
				   const SystShifts& /*syst*/) const
  {
    return util::sqr((fVar->GetValue(osc)-fMean)/fSigma);
  }

}
