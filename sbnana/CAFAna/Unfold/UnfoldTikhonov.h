#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"

namespace ana
{
  /// \brief Unfolding using Tikhonov regularization (penalizing true spectra
  /// with large second derivatives)
  ///
  /// Currently requires 1D reco and true axes, but that is not a fundamental
  /// restriction of the method and could be lifted in future versions.
  ///
  /// Large \a regStrength corresponds to more regularization (biased towards
  /// smooth results), small strength suffers more from statistical
  /// fluctuations.
  Spectrum UnfoldTikhonov(const Spectrum& reco,
                          const ReweightableSpectrum& recoVsTrue,
                          double regStrength);
}
