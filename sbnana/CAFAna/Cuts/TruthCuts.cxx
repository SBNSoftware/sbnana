#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  const Cut kHasMatchedNu([](const caf::SRSliceProxy* slc)
                          {
                            return slc->truth.index >= 0;
                          });

  /// \brief Is this a Neutral %Current event?
  ///
  /// In this case the selection function is simple enough that we can include
  /// it inline as a lambda function.
  const Cut kIsNC([](const caf::SRSliceProxy* slc)
                  {
                    return kHasMatchedNu(slc) && slc->truth.isnc;
                  });

  const Cut kIsCC([](const caf::SRSliceProxy* slc)
                  {
                    return kHasMatchedNu(slc) && slc->truth.iscc;
                  });

  const Cut kIsAntiNu([](const caf::SRSliceProxy* slc)
                      {
                        return kHasMatchedNu(slc) && slc->truth.pdg < 0;
                      });

  //----------------------------------------------------------------------
  bool CCFlavSel::operator()(const caf::SRSliceProxy* slc) const
  {
    return kHasMatchedNu(slc) &&
      slc->truth.iscc &&
      abs(slc->truth.initpdg) == fPdgOrig &&
      abs(slc->truth.pdg) == fPdg;
  }

  //----------------------------------------------------------------------
  bool NCFlavOrig::operator()(const caf::SRSliceProxy* slc) const
  {
    return kHasMatchedNu(slc) && slc->truth.isnc && abs(slc->truth.initpdg) == fPdgOrig;
  }
}
