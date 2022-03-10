#pragma once

#include "cafanacore/MultiVar.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  typedef _MultiVar<caf::SRSliceProxy> MultiVar;
  typedef _MultiVar<caf::SRSpillProxy> SpillMultiVar;
} // namespace
