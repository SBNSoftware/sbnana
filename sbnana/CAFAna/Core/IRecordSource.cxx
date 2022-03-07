#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  //----------------------------------------------------------------------
  class ShiftedSliceEnsembleSource:
    public beta::PassthroughExposure<ISliceSink,
                                     ISliceEnsembleSource>
  {
  public:
    ShiftedSliceEnsembleSource(ISliceSource& src,
                               const std::vector<SystShifts>& shifts,
                               int multiverseId);

    virtual void HandleRecord(const caf::SRSliceProxy* slc, double weight) override;

  protected:
    std::vector<SystShifts> fShifts;
    int fMultiverseId;
  };

  //----------------------------------------------------------------------
  ShiftedSliceEnsembleSource::
  ShiftedSliceEnsembleSource(ISliceSource& src,
                             const std::vector<SystShifts>& shifts,
                             int multiverseId)
    : fShifts(shifts), fMultiverseId(multiverseId)
  {
    src.Register(this);
  }

  //----------------------------------------------------------------------
  void ShiftedSliceEnsembleSource::HandleRecord(const caf::SRSliceProxy* slc,
                                                double weight)
  {
    if(weight == 0) return;

    std::vector<double> weights(fShifts.size(), weight);

    bool anyShifted = false;

    for(unsigned int univIdx = 0; univIdx < fShifts.size(); ++univIdx){
      // Need to provide a clean slate for each new set of systematic
      // shifts to work from. Copying the whole StandardRecord is pretty
      // expensive, so modify it in place and revert it afterwards.

      caf::SRProxySystController::BeginTransaction();

      bool shifted = false;

      // Can special-case nominal to not pay cost of Shift()
      if(!fShifts[univIdx].IsNominal()){
        // const_cast is naughty. I hope the fact that we put the record back
        // afterwards absolves most sins.
        fShifts[univIdx].Shift(const_cast<caf::SRSliceProxy*>(slc), weights[univIdx]);
        // If there were only weighting systs applied then the cached
        // nominal values are still valid.
        shifted = caf::SRProxySystController::AnyShifted();
      }

      // Slice was shifted or we are already in the slow path, so we have to
      // handle this individually
      if((shifted || anyShifted) && weights[univIdx] != 0)
        for(ISliceEnsembleSink* sink: fSinks)
          sink->HandleSingleRecord(slc, weights[univIdx], univIdx);

      // Return StandardRecord to its unshifted form ready for the next
      // histogram.
      caf::SRProxySystController::Rollback();

      // We entered the slow path for the first time for this universe, have to
      // catch up with all the other universes we were hoping to be able to
      // handle in the fast path.
      if(shifted && !anyShifted){
        anyShifted = true;

        for(unsigned int prevIdx = 0; prevIdx < univIdx; ++prevIdx){
          if(weights[prevIdx] == 0) continue;

          for(ISliceEnsembleSink* sink: fSinks){
            sink->HandleSingleRecord(slc, weights[prevIdx], prevIdx);
          }
        } // end for prevIdx
      } // end if shifted for the first time
    } // end for univIdx

    // Fast path in case none of the records got rewritten, can treat as an
    // ensemble with weights.
    if(!anyShifted){
      for(ISliceEnsembleSink* sink: fSinks){
        sink->HandleEnsemble(slc, weights);
      }
    }
  }

  //----------------------------------------------------------------------
  ISliceEnsembleSource& ISliceSource::
  Ensemble(const std::vector<SystShifts>& shifts, int multiverseId)
  {
    return fEnsembleSources.template Get<ShiftedSliceEnsembleSource>(multiverseId, *this, shifts, multiverseId);
  }

  //----------------------------------------------------------------------
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
