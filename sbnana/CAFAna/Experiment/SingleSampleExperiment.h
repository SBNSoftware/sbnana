#pragma once

#include "sbnana/CAFAna/Experiment/IExperiment.h"
#include "sbnana/CAFAna/Prediction/IPrediction.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

namespace ana
{
  /// Compare a single data spectrum to the MC + cosmics expectation
  class SingleSampleExperiment: public IExperiment
  {
  public:
    /// \param pred            Source of oscillated MC beam predictions
    /// \param data            Data spectrum to compare to
    /// \param cosmicInTime    In-time cosmic ray background component
    /// \param cosmicOutOfTime Out-of-time cosmic ray background component
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data,
                           const Spectrum& cosmicInTime,
                           const Spectrum& cosmicOutOfTime);

    /// In MC studies you might not want to bother with cosmics
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data)
      : fMC(pred), fData(data), fCosmicInTime(0, {}, {}, 0, 0), fCosmicOutOfTime(0, {}, {}, 0, 0), fMask(0)
    {
    }

    virtual ~SingleSampleExperiment();

    virtual double ChiSq(osc::IOscCalcAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<SingleSampleExperiment> LoadFrom(TDirectory* dir);

    // Didn't make provisions for copying fCosmic or fMC
    SingleSampleExperiment(const SingleSampleExperiment&) = delete;
    SingleSampleExperiment& operator=(const SingleSampleExperiment&) = delete;

    // need to explicitly declare move constructor since copy constructor is deleted
    SingleSampleExperiment(SingleSampleExperiment&& s)
      : fMC(s.fMC),
        fData(std::move(s.fData)),
        fCosmicInTime(std::move(s.fCosmicInTime)),
        fCosmicOutOfTime(std::move(s.fCosmicOutOfTime))
    {
      s.fMC = nullptr;
    };

    void SetMaskHist(double xmin=0, double xmax=-1, 
		     double ymin=0, double ymax=-1);

  protected:
    const IPrediction* fMC;
    Spectrum fData;
    Spectrum fCosmicInTime, fCosmicOutOfTime;
    TH1* fMask;
  };
}
