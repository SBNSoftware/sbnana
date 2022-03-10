#pragma once

#include "cafanacore/HistAxis.h"

#include "sbnana/CAFAna/Core/Var.h" // TODO do we want our own FwdDeclare.h?

namespace ana
{
  typedef _HistAxis<Var> HistAxis;
  typedef _HistAxis<SpillVar> SpillHistAxis;

  typedef _HistAxis<TrackVar> TrackHistAxis;
  typedef _HistAxis<ShowerVar> ShowerHistAxis;
  typedef _HistAxis<StubVar> StubHistAxis;
}
