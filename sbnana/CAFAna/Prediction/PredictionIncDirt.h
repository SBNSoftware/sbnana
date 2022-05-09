#pragma once

#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"

//#include "sbnana/CAFAna/Prediction/PredictionGenerator.h"

namespace ana
{
  class Loaders;

  /// Prediction summing detector and dirt components
  class PredictionIncDirt: public IPrediction
  {
  public:
    PredictionIncDirt(ISliceSource& srcNonswap,
                      ISliceSource& srcNue,
                      ISliceSource& srcNuTau,
                      ISliceSource& srcIntrinsic,
                      ISliceSource& srcDirt,
                      const HistAxis& axis);

    PredictionIncDirt(SliceSources& srcs,
                      ISliceSource& srcDirt,
                      const HistAxis& axis);

    virtual ~PredictionIncDirt();

    static std::unique_ptr<PredictionIncDirt> LoadFrom(TDirectory* dir, const std::string& name);
    virtual void SaveTo(TDirectory* dir, const std::string& name) const override;

    Spectrum PredictDet(osc::IOscCalc* calc) const
    {
      return fDet.Predict(calc);
    }

    Spectrum PredictDirt(osc::IOscCalc* calc) const
    {
      return fDirt.Predict(calc);
    }

    Spectrum PredictComponentDet(osc::IOscCalc* calc,
                                 Flavors::Flavors_t flav,
                                 Current::Current_t curr,
                                 Sign::Sign_t sign) const
    {
      return fDet.PredictComponent(calc, flav, curr, sign);
    }

    Spectrum PredictComponentDirt(osc::IOscCalc* calc,
                                  Flavors::Flavors_t flav,
                                  Current::Current_t curr,
                                  Sign::Sign_t sign) const
    {
      return fDirt.PredictComponent(calc, flav, curr, sign);
    }
    
    virtual Spectrum Predict(osc::IOscCalc* calc) const override
    {
      return PredictDet(calc) + PredictDirt(calc);
    }

    virtual Spectrum PredictComponent(osc::IOscCalc* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override
    {
      return (PredictComponentDet (calc, flav, curr, sign) +
              PredictComponentDirt(calc, flav, curr, sign));
    }

  protected:
    PredictionIncDirt(std::unique_ptr<PredictionNoExtrap>&& det,
                      std::unique_ptr<PredictionNoExtrap>&& dirt)
      : fDet(*det), fDirt(*dirt)
    {
    }

    PredictionNoExtrap fDet, fDirt;
  };

  // TODO how best to write a PredictionGenerator for this?
}
