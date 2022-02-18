#pragma once

#include "CAFAnaCore/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/IRecordSink.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  using ISpillSource = beta::_IRecordSource<caf::SRSpillProxy>;
  using ISliceSource = beta::_IRecordSource<caf::SRSliceProxy>;

  using ISpillEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSpillProxy>;
  using ISliceEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSliceProxy>;



  class SliceAdaptor : public beta::PassthroughExposure<ISpillSink, ISliceSource>
  {
  public:
    SliceAdaptor(ISpillSource& src);

    virtual void HandleRecord(const caf::SRSpillProxy* spill, double weight) override;

    // Will need an EnsembleSliceAdaptor if we ever have syst weights that
    // apply at the spill level.

    //    virtual void HandleEnsemble(const caf::SRSpillProxy* spill, const std::vector<double>& weights, int multiverseId) override;
  };


  // Spill sources are also slice sources (they just loop over the slices)
  template<> class beta::_IRecordSource<caf::SRSpillProxy> : public beta::_IRecordSourceDefaultImpl<caf::SRSpillProxy>
  {
  public:
    ISliceSource& Slices() {return fSlices;}

  protected:
    _IRecordSource() : fSlices(*this) {}

    SliceAdaptor fSlices;
  };
}
