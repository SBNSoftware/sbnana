#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"

namespace ana
{
  Spectrum UnfoldTikhonov(const Spectrum& reco,
                          const ReweightableSpectrum& recoVsTrue,
                          double regStrength);
}
