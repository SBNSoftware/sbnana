#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"

namespace ana
{
  Spectrum UnfoldSVD(const Spectrum& reco,
                     const ReweightableSpectrum& recoVsTrue,
                     unsigned int reg);
}
