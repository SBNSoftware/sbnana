#pragma once

#include <memory>

#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Core/Loaders.h"

namespace ana
{

  class Loaders;
  class IPrediction;

  /** \brief Given loaders and an MC shift, Generate() generates an IPrediction

      All other arguments needed to construct the prediction are passed to the
      IPredictionGenerator constructor, and are identical for all generated
      predictions. For standalone use or use with PredictionInterp. **/
  class IPredictionGenerator
  {
    public:
    virtual std::unique_ptr<IPrediction> Generate(
    						  Loaders& loaders, const SystShifts& shiftMC = kNoShift ) const = 0;

  };

}
