#pragma once

#include "sbnana/CAFAna/Experiment/IExperiment.h"
#include "sbnana/CAFAna/Prediction/IPrediction.h"

#include "CAFAna/Core/FwdDeclare.h"

namespace ana
{
  /// Compare a single data spectrum to the MC + cosmics expectation
  class SingleSampleExperiment: public IExperiment
  {
  public:
    /// \param pred   Source of oscillated MC beam predictions
    /// \param data   Data spectrum to compare to
    /// \param cosmic In-time cosmic ray background component
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data,
                           const Spectrum& cosmic);

    /// In MC studies you might not want to bother with cosmics
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data)
      : fMC(pred), fData(data), fCosmic(Spectrum::Uninitialized())
    {
    }

    virtual ~SingleSampleExperiment();

    virtual double ChiSq(osc::IOscCalcAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<SingleSampleExperiment> LoadFrom(TDirectory* dir);

    // Didn't make provisions for copying fMC
    SingleSampleExperiment(const SingleSampleExperiment&) = delete;
    SingleSampleExperiment& operator=(const SingleSampleExperiment&) = delete;

    // need to explicitly declare move constructor since copy constructor is deleted
    SingleSampleExperiment(SingleSampleExperiment&& s)
      : fMC(s.fMC), fData(std::move(s.fData)), fCosmic(std::move(s.fCosmic))
    {
      s.fMC = nullptr;
    };

    void SetMaskHist(double xmin=0, double xmax=-1, 
		     double ymin=0, double ymax=-1);

  protected:
    virtual void ApplyMask(Eigen::ArrayXd& a, Eigen::ArrayXd& b) const;

    const IPrediction* fMC;
    Spectrum fData;
    Spectrum fCosmic;

    Eigen::ArrayXd fMaskA;
  };
}
