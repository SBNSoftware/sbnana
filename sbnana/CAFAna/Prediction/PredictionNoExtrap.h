#pragma once

#include "sbnana/CAFAna/Prediction/PredictionExtrap.h"

#include "sbnana/CAFAna/Prediction/PredictionGenerator.h"

namespace ana
{
  class Loaders;

  /// Prediction that just uses one detector's MC, with no extrapolation
  class PredictionNoExtrap: public PredictionExtrap
  {
  public:
    PredictionNoExtrap(PredictionExtrap* pred);

    PredictionNoExtrap(SpectrumLoaderBase& loaderNonswap,
                       SpectrumLoaderBase& loaderNue,
                       SpectrumLoaderBase& loaderNuTau,
                       SpectrumLoaderBase& loaderIntrinsic,
                       const std::string& label,
                       const Binning& bins,
                       const Var& var,
                       const SpillCut& spillcut,
                       const Cut& cut,
                       const SystShifts& shift = kNoShift,
                       const Var& wei = kUnweighted);

    PredictionNoExtrap(SpectrumLoaderBase& loaderNonswap,
                       SpectrumLoaderBase& loaderNue,
                       SpectrumLoaderBase& loaderNuTau,
                       SpectrumLoaderBase& loaderIntrinsic,
		       const HistAxis& axis,
                       const SpillCut& spillcut,
		       const Cut& cut,
                       const SystShifts& shift = kNoShift,
                       const Var& wei = kUnweighted);

    PredictionNoExtrap(Loaders& loaders,
                       const std::string& label,
                       const Binning& bins,
                       const Var& var,
                       const SpillCut& spillcut,
                       const Cut& cut,
                       const SystShifts& shift = kNoShift,
                       const Var& wei = kUnweighted);

    PredictionNoExtrap(Loaders& loaders,
                       const HistAxis& axis,
                       const SpillCut& spillcut,
                       const Cut& cut,
                       const SystShifts& shift = kNoShift,
                       const Var& wei = kUnweighted);

    virtual ~PredictionNoExtrap();

    static std::unique_ptr<PredictionNoExtrap> LoadFrom(TDirectory* dir);
    virtual void SaveTo(TDirectory* dir) const override;

  };

  class NoExtrapPredictionGenerator: public IPredictionGenerator
  {
  public:
    NoExtrapPredictionGenerator(HistAxis axis, SpillCut spillcut, Cut cut, Var wei = kUnweighted)
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
    Var fWei;
  };
}
