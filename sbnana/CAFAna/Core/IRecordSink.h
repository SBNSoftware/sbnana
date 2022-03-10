#pragma once

#include "cafanacore/IRecordSink.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  using ISpillSink = beta::_IRecordSink<caf::SRSpillProxy>;
  using ISliceSink = beta::_IRecordSink<caf::SRSliceProxy>;

  using ISpillEnsembleSink = beta::_IRecordEnsembleSink<caf::SRSpillProxy>;
  using ISliceEnsembleSink = beta::_IRecordEnsembleSink<caf::SRSliceProxy>;
}
