#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"

namespace ana
{
  // \brief Iterative unfolding, also known as "Bayesian" or "d'Agostini"
  //
  // Smaller values of \a nIterations are more regularized (more biased to the
  // MC distribution), large values suffer more from fluctuations.
  Spectrum UnfoldIterative(const Spectrum& reco,
                           const ReweightableSpectrum& recoVsTrue,
                           unsigned int nIterations);
}
