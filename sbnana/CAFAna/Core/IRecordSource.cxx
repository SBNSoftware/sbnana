#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  //----------------------------------------------------------------------
  SliceAdaptor::SliceAdaptor(ISpillSource& src)
  {
    src.Register(this);
  }

  //----------------------------------------------------------------------
  void SliceAdaptor::HandleRecord(const caf::SRSpillProxy* spill,
                                  double weight)
  {
    for(const caf::SRSliceProxy& slc: spill->slc)
      for(auto& sink: fSinks)
        sink->HandleRecord(&slc, weight);
  }

  /*
  //----------------------------------------------------------------------
  void SliceAdaptor::HandleEnsemble(const caf::SRSpillProxy* spill,
                                    const std::vector<double>& weights,
                                    int multiverseId)
  {
    for(const caf::SRSliceProxy& slc: spill->slc)
      for(auto& sink: fSinks)
        sink->HandleEnsemble(&slc, weights, multiverseId);
  }
  */
}
