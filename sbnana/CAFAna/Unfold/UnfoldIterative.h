#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"

namespace ana
{
  Spectrum UnfoldIterative(const Spectrum& reco,
                           const ReweightableSpectrum& recoVsTrue,
                           unsigned int nIterations);
}
