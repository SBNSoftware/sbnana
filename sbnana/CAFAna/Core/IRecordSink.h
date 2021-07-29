#include "CAFAnaCore/CAFAna/Core/IRecordSink.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  using ISpillSink = beta::_IRecordSink<caf::SRSpillProxy>;
  using ISliceSink = beta::_IRecordSink<caf::SRSliceProxy>;
}
