#pragma once

#include "CAFAnaCore/CAFAna/Core/Weight.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  typedef _Weight<caf::SRSliceProxy> Weight;

  /// \brief Equivalent of \ref Weight acting on \ref caf::SRSpill
  typedef _Weight<caf::SRSpillProxy> SpillWeight;

  /// The simplest possible Weight, always 1. Used as a default weight.
  const Weight kUnweighted = Unweighted<caf::SRSliceProxy>();

  const SpillWeight kSpillUnweighted = Unweighted<caf::SRSpillProxy>();
} // namespace
