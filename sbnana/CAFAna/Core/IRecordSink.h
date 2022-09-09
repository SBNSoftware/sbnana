#pragma once

#include "cafanacore/IRecordSink.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  using ISpillSink = _IRecordSink<caf::SRSpillProxy>;
  using ISliceSink = _IRecordSink<caf::SRSliceProxy>;

  using ISpillEnsembleSink = _IRecordEnsembleSink<caf::SRSpillProxy>;
  using ISliceEnsembleSink = _IRecordEnsembleSink<caf::SRSliceProxy>;
}
