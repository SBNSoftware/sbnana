#pragma once

#include "sbnana/CAFAna/Experiment/IExperiment.h"
#include "CAFAna/Core/Spectrum.h"


namespace ana
{
  class IPrediction;

  /// Compare a data spectrum to MC expectation using only the event count
  class CountingExperiment: public IExperiment
  {
  public:
    CountingExperiment(const IPrediction* p, const Spectrum& d) : fMC(p), fData(d) {}
    ~CountingExperiment();
    virtual double ChiSq(osc::IOscCalcAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<CountingExperiment> LoadFrom(TDirectory* dir);
  protected:
    const IPrediction* fMC;
    Spectrum fData;
  };
}
