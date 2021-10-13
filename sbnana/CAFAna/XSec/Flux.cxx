#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  //----------------------------------------------------------------------
  bool IsCCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg)
  {
    // It doesn't matter exactly what process we choose, so long as it
    // corresponds to one of the modes GENIE considers as a top-level
    // cross-section. Here we go for CC QE on Argon.
    return nu->pdg == pdg && nu->iscc &&
      nu->genie_mode == caf::kQE && nu->genie_inttype == caf::kCCQE &&
      nu->targetPDG == 1000180400 /* Argon 40 */ &&
      !nu->ischarm;
  }

  //----------------------------------------------------------------------
  bool IsNCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg)
  {
    // This process works for very low neutrino energy where the muon mass
    // becomes relevant.
    return nu->pdg == pdg && !nu->iscc &&
      nu->genie_mode == caf::kQE && nu->genie_inttype == caf::kNCQE &&
      nu->targetPDG == 1000180400 /* Argon 40 */ && 
      !nu->ischarm &&
      nu->hitnuc == 2112;
  }

}
