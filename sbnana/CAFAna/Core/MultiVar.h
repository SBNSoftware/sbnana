#pragma once

#include "CAFAnaCore/CAFAna/Core/MultiVar.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  typedef _MultiVar<caf::SRSliceProxy> MultiVar;
  typedef _MultiVar<caf::SRSpillProxy> SpillMultiVar;
} // namespace
