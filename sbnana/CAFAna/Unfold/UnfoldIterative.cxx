#include "sbnana/CAFAna/Unfold/UnfoldIterative.h"

#include <cassert>

namespace ana
{
  //----------------------------------------------------------------------
  Spectrum UnfoldIterative(const Spectrum& reco,
                           const ReweightableSpectrum& recoVsTrue,
                           unsigned int nIterations)
  {
    assert(nIterations > 0);

    // Working copy
    ReweightableSpectrum rt = recoVsTrue;

    rt.ReweightToRecoSpectrum(reco);

    // First result after one iteration
    Spectrum truth = rt.WeightingVariable();

    for(unsigned int i = 2; i <= nIterations; ++i){
      rt = recoVsTrue; // put it back how it started
      rt.ReweightToTrueSpectrum(truth); // project back the other way
      rt.ReweightToRecoSpectrum(reco); // and then forward again
      truth = rt.WeightingVariable();
    }

    return truth.FakeData(reco.POT());
  }
}
