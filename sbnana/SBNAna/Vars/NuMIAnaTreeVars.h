#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana
{

  /// Which type of selection is this: Signal (1), Chg Pi sideband (2), Pi0 sideband (3), neither (0)
  extern const Var kNuMICutType;

  /// Which type of event is this: Signal (1), Other NuCC (2), NC (3), Cosmic/Unmatched (4)
  extern const Var kNuMISliceCategory;

  /// Which true neutrino pdg is this slice matching to, if any?
  extern const Var kNuMITruePDG;

  /// Is it FHC or RHC -- NB: for now, this is only returning 1 (== FHC)
  extern const Var kNuMIIsFHC;

  /// Interaction mode from GENIE
  extern const Var kNuMITrueMode;

  /// True target nu interacts with
  extern const Var kNuMITrueTarget;

  /// Is the truth matched interaction charged-current?
  extern const Var kNuMITrueIsCC;

  /// True neutrino E
  extern const Var kNuMITrueENu;

}
