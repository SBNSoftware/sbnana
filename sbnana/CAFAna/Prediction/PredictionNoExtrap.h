#pragma once

#include "sbnana/CAFAna/Prediction/PredictionExtrap.h"

#include "sbnana/CAFAna/Prediction/PredictionGenerator.h"

#include "sbnana/CAFAna/Extrap/TrivialExtrap.h"

namespace ana
{
  class Loaders;

  /// Prediction that just uses one detector's MC, with no extrapolation
  class PredictionNoExtrap: public PredictionExtrap
  {
  public:
    PredictionNoExtrap(ISliceSource& srcNonswap,
                       ISliceSource& srcNue,
                       ISliceSource& srcNuTau,
                       ISliceSource& srcIntrinsic,
                       const HistAxis& axis);

    PredictionNoExtrap(SliceSources& srcs, const HistAxis& axis);

    virtual ~PredictionNoExtrap();

    static std::unique_ptr<PredictionNoExtrap> LoadFrom(TDirectory* dir, const std::string& name);
    virtual void SaveTo(TDirectory* dir, const std::string& name) const override;

  protected:
    PredictionNoExtrap(TrivialExtrap* extrap) : PredictionExtrap(extrap) {}
  };

  /* TODO think about generators
  class NoExtrapPredictionGenerator: public IPredictionGenerator
  {
  public:
    NoExtrapPredictionGenerator(HistAxis axis, SpillCut spillcut, Cut cut, Weight wei = kUnweighted)
      : fAxis(axis), fSpillCut(spillcut), fCut(cut), fWei(wei)
    {
    }

    virtual std::unique_ptr<IPrediction>
    Generate(Loaders& loaders, const SystShifts& shiftMC = kNoShift) const override
    {
      return std::unique_ptr<IPrediction>(new PredictionNoExtrap(loaders, fAxis, fSpillCut, fCut, shiftMC, fWei));
    }

  protected:
    HistAxis fAxis;
    SpillCut fSpillCut;
    Cut fCut;
    Weight fWei;
  };
  */
}
