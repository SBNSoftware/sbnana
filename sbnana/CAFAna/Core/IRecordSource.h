#include "CAFAnaCore/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  using ISpillSource = beta::_IRecordSource<caf::SRSpillProxy>;
  using ISliceSource = beta::_IRecordSource<caf::SRSliceProxy>;
}
