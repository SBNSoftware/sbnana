#pragma once

#include "cafanacore/Cut.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  /// \brief Representation of a cut (selection) to be applied to a \ref
  /// caf::StandardRecord object
  ///
  /// A Cut consists of a function, taking a StandardRecord and returning a
  /// boolean indicating if that event passes the cut.
  ///
  /// Cut objects may be combined with the standard boolean operations && ||
  /// and !
  typedef _Cut<caf::SRSliceProxy> SliceCut;
  typedef _Cut<caf::SRSliceProxy> Cut;

  /// \brief Equivalent of \ref Cut acting on \ref caf::SRSpill. For use in
  /// spill-by-spill data quality cuts
  typedef _Cut<caf::SRSpillProxy> SpillCut;

  /// The simplest possible cut: pass everything, used as a default
  const Cut kNoCut(NoCut<caf::SRSliceProxy>{});

  /// The simplest possible cut: pass everything, used as a default
  const SpillCut kNoSpillCut(NoCut<caf::SRSpillProxy>{});


  typedef _Cut<caf::SRTrueInteractionProxy> NuTruthCut;
  typedef _Cut<caf::SRTrackProxy> TrackCut;
  typedef _Cut<caf::SRShowerProxy> ShowerCut;
  typedef _Cut<caf::SRStubProxy> StubCut;

  const NuTruthCut  kNoNuTruthCut(NoCut<caf::SRTrueInteractionProxy>{});
  const TrackCut kNoTrackCut(NoCut<caf::SRTrackProxy>{});
  const ShowerCut kNoShowerCut(NoCut<caf::SRShowerProxy>{});
  const StubCut kNoStubCut(NoCut<caf::SRStubProxy>{});
} // namespace
