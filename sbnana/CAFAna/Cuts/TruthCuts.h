#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  extern const Cut kHasMatchedNu;

  /// \brief Is this a Neutral %Current event?
  ///
  /// In this case the selection function is simple enough that we can include
  /// it inline as a lambda function.
  extern const Cut kIsNC;

  extern const Cut kIsCC;

  //----------------------------------------------------------------------
  /// Helper for defining true CC event cuts
  class CCFlavSel
  {
  public:
    CCFlavSel(int pdg, int pdgorig) : fPdg(pdg), fPdgOrig(pdgorig)
    {
    }

    bool operator()(const caf::SRSliceProxy* slc) const;
  protected:
    int fPdg, fPdgOrig;
  };

  /// Helper for defining true CC event cuts
  class NCFlavOrig
  {
  public:
    NCFlavOrig(int pdgorig) : fPdgOrig(pdgorig)
    {
    }

    bool operator()(const caf::SRSliceProxy* slc) const;
  protected:
    int fPdgOrig;
  };

  // Finally, the function argument to the Cut constructor can be a "functor"
  // object (one with operator()). This allows similar logic but with different
  // constants to be easily duplicated.

  /// Select CC \f$ \nu_\mu\to\nu_e \f$
  const Cut kIsNueApp (CCFlavSel(12, 14));
  /// Select CC \f$ \nu_\mu\to\nu_\mu \f$
  const Cut kIsNumuCC (CCFlavSel(14, 14));
  /// Select CC \f$ \nu_e\to\nu_e \f$
  const Cut kIsBeamNue(CCFlavSel(12, 12));
  /// Select CC \f$ \nu_e\to\nu_\mu \f$
  const Cut kIsNumuApp(CCFlavSel(14, 12));
  /// Select CC \f$ \nu_\mu\to\nu_\tau \f$
  const Cut kIsTauFromMu(CCFlavSel(16, 14));
  /// Select CC \f$ \nu_e\to\nu_\tau \f$
  const Cut kIsTauFromE(CCFlavSel(16, 12));


  // NC from Numu
  const Cut kIsNCFromNumu(NCFlavOrig(14));
  // NC from Nue
  const Cut kIsNCFromNue(NCFlavOrig(12));

  /// Is this truly an antineutrino?
  extern const Cut kIsAntiNu;
}
