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
  TrackAdaptor::TrackAdaptor(ISliceRecoBranchSource& src)
  {
    src.Register(this);
  }

  //----------------------------------------------------------------------
  ShowerAdaptor::ShowerAdaptor(ISliceRecoBranchSource& src)
  {
    src.Register(this);
  }

  //----------------------------------------------------------------------
  StubAdaptor::StubAdaptor(ISliceRecoBranchSource& src)
  {
    src.Register(this);
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  void SliceAdaptor::HandleRecord(const caf::SRSpillProxy* spill,
                                  double weight)
  {
    for(const caf::SRSliceProxy& slc: spill->slc)
      for(auto& sink: fSinks)
        sink->HandleRecord(&slc, weight);
  }

  //----------------------------------------------------------------------
  void TrackAdaptor::HandleRecord(const caf::SRSliceRecoBranchProxy* reco,
                                  double weight)
  {
    for(const caf::SRTrackProxy& trk: reco->trk)
      for(auto& sink: fSinks)
        sink->HandleRecord(&trk, weight);
  }

  //----------------------------------------------------------------------
  void ShowerAdaptor::HandleRecord(const caf::SRSliceRecoBranchProxy* reco,
                                   double weight)
  {
    for(const caf::SRShowerProxy& shw: reco->shw)
      for(auto& sink: fSinks)
        sink->HandleRecord(&shw, weight);
  }

  //----------------------------------------------------------------------
  void StubAdaptor::HandleRecord(const caf::SRSliceRecoBranchProxy* reco,
                                 double weight)
  {
    for(const caf::SRStubProxy& stub: reco->stub)
      for(auto& sink: fSinks)
        sink->HandleRecord(&stub, weight);
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  void TrackEnsembleAdaptor::
  HandleSingleRecord(const caf::SRSliceRecoBranchProxy* reco,
                     double weight,
                     int universeIdx)
  {
    for(const caf::SRTrackProxy& trk: reco->trk)
      for(auto& sink: fSinks)
        sink->HandleSingleRecord(&trk, weight, universeIdx);
  }

  //----------------------------------------------------------------------
  void ShowerEnsembleAdaptor::
  HandleSingleRecord(const caf::SRSliceRecoBranchProxy* reco,
                     double weight,
                     int universeIdx)
  {
    for(const caf::SRShowerProxy& shw: reco->shw)
      for(auto& sink: fSinks)
        sink->HandleSingleRecord(&shw, weight, universeIdx);
  }

  //----------------------------------------------------------------------
  void StubEnsembleAdaptor::
  HandleSingleRecord(const caf::SRSliceRecoBranchProxy* reco,
                     double weight,
                     int universeIdx)
  {
    for(const caf::SRStubProxy& stub: reco->stub)
      for(auto& sink: fSinks)
        sink->HandleSingleRecord(&stub, weight, universeIdx);
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  void TrackEnsembleAdaptor::
  HandleEnsemble(const caf::SRSliceRecoBranchProxy* reco,
                 const std::vector<double>& weights)
  {
    for(const caf::SRTrackProxy& trk: reco->trk)
      for(auto& sink: fSinks)
        sink->HandleEnsemble(&trk, weights);
  }

  //----------------------------------------------------------------------
  void ShowerEnsembleAdaptor::
  HandleEnsemble(const caf::SRSliceRecoBranchProxy* reco,
                 const std::vector<double>& weights)
  {
    for(const caf::SRShowerProxy& shw: reco->shw)
      for(auto& sink: fSinks)
        sink->HandleEnsemble(&shw, weights);
  }

  //----------------------------------------------------------------------
  void StubEnsembleAdaptor::
  HandleEnsemble(const caf::SRSliceRecoBranchProxy* reco,
                 const std::vector<double>& weights)
  {
    for(const caf::SRStubProxy& stub: reco->stub)
      for(auto& sink: fSinks)
        sink->HandleEnsemble(&stub, weights);
  }

}
