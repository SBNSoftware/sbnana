#include "sbnana/CAFAna/Prediction/IPrediction.h"
#include "sbnana/CAFAna/Prediction/PredictionGenerator.h"

namespace ana
{
  class Loaders;

  /// Prediction that wraps a simple Spectrum
  class PredictionNoOsc: public IPrediction
  {
  public:
    PredictionNoOsc(ISliceSource& src, const HistAxis& axis);

    static std::unique_ptr<PredictionNoOsc> LoadFrom(TDirectory* dir);
    virtual void SaveTo(TDirectory* dir) const override;

    virtual Spectrum Predict(osc::IOscCalc* /*calc*/) const override
    {
      return fSpectrum;
    }

    virtual Spectrum PredictComponent(osc::IOscCalc* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override;

  protected:
    PredictionNoOsc(const Spectrum& s,
                    const Spectrum& sNC,
                    const Spectrum& sNumu, const Spectrum& sNumubar,
                    const Spectrum& sNue, const Spectrum& sNuebar)
      : fSpectrum(s),
        fSpectrumNC(sNC),
        fSpectrumNumu(sNumu), fSpectrumNumubar(sNumubar),
        fSpectrumNue(sNue), fSpectrumNuebar(sNuebar)
    {
    }

    Spectrum fSpectrum;

    Spectrum fSpectrumNC;
    Spectrum fSpectrumNumu;
    Spectrum fSpectrumNumubar;
    Spectrum fSpectrumNue;
    Spectrum fSpectrumNuebar;
  };


  class NoOscPredictionGenerator: public IPredictionGenerator
  {
  public:
    NoOscPredictionGenerator(SpectrumLoaderBase& loader,
                             HistAxis axis,
                             SpillCut spillcut,
                             Cut cut,
			     Weight wei = kUnweighted)
      : fLoader(loader), fAxis(axis), fSpillCut(spillcut), fCut(cut), fWei(wei)
    {
    }

    virtual std::unique_ptr<IPrediction>
    Generate(Loaders& loaders, const SystShifts& shiftMC = kNoShift) const override
    {
      abort();
      // TODO TODO TODO
      /*
      return std::unique_ptr<IPrediction>(new PredictionNoOsc(fLoader,
                                                              fAxis,
                                                              fSpillCut,
                                                              fCut,
                                                              shiftMC, fWei));
      */
    }
  protected:
    SpectrumLoaderBase& fLoader;
    HistAxis fAxis;
    SpillCut fSpillCut;
    Cut fCut;
    Weight fWei;
  };

}
