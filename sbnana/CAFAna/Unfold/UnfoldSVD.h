#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"

namespace ana
{
  /// \brief Singular value decomposition unfolding using ROOT's TSVDDecomp
  ///
  /// SVD decomposition requires the reco and true spectra to have the same
  /// number of bins (recoVsTrue must be square).
  ///
  /// Lower values of \a reg are more regularized
  Spectrum UnfoldSVD(const Spectrum& reco,
                     const ReweightableSpectrum& recoVsTrue,
                     unsigned int reg);
}
