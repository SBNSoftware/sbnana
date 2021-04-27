#include "sbnana/CAFAna/Prediction/PredictionGenerator.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"

namespace ana
{

  //--------------------------------------------------------------------------

  NoExtrapGenerator::NoExtrapGenerator(
    const HistAxis axis,
    const SpillCut spillcut,
    const Cut cut,
    const Var wei
  ) : fAxis(axis), fSpillCut(spillcut), fCut(cut), fWei(wei) {}

  std::unique_ptr<IPrediction> NoExtrapGenerator::Generate(
							   Loaders& loaders,
							   const SystShifts& shiftMC
							   ) const {
    return std::unique_ptr<IPrediction>( new PredictionNoExtrap(
								loaders, fAxis, fSpillCut, fCut, shiftMC, fWei ) );
  }

}
