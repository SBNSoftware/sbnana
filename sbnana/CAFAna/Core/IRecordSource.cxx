#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/SystShifts.h"

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

  template<class FromT, class ToT> VectorAdaptor<FromT, ToT>::
  VectorAdaptor(beta::_IRecordSource<caf::Proxy<FromT>>& src,
                Func_t vecGetter)
    : fVecGetter(vecGetter)
  {
    src.Register(this);
  }

  //----------------------------------------------------------------------
  template<class FromT, class ToT> void VectorAdaptor<FromT, ToT>::
  HandleRecord(const caf::Proxy<FromT>* rec, double weight)
  {
    for(const caf::Proxy<ToT>& to: fVecGetter(rec))
      for(auto& sink: beta::_IRecordSource<caf::Proxy<ToT>>::fSinks)
        sink->HandleRecord(&to, weight);
  }

  //----------------------------------------------------------------------
  template<class FromT, class ToT> EnsembleVectorAdaptor<FromT, ToT>::
  EnsembleVectorAdaptor(beta::_IRecordEnsembleSource<caf::Proxy<FromT>>& src,
                        Func_t vecGetter)
    : fVecGetter(vecGetter)
  {
    src.Register(this);
  }

  //----------------------------------------------------------------------
  template<class FromT, class ToT> void EnsembleVectorAdaptor<FromT, ToT>::
  HandleSingleRecord(const caf::Proxy<FromT>* rec,
                     double weight,
                     int universeId)
  {
    for(const caf::Proxy<ToT>& to: fVecGetter(rec))
      for(auto& sink: beta::_IRecordEnsembleSource<caf::Proxy<ToT>>::fSinks)
        sink->HandleSingleRecord(&to, weight, universeId);
  }

  //----------------------------------------------------------------------
  template<class FromT, class ToT> void EnsembleVectorAdaptor<FromT, ToT>::
  HandleEnsemble(const caf::Proxy<FromT>* rec,
                 const std::vector<double>& weights)
  {
    for(const caf::Proxy<ToT>& to: fVecGetter(rec))
      for(auto& sink: beta::_IRecordEnsembleSource<caf::Proxy<ToT>>::fSinks)
        sink->HandleEnsemble(&to, weights);
  }

  //----------------------------------------------------------------------
  const caf::Proxy<std::vector<caf::SRSlice>>&
  GetSlices(const caf::SRSpillProxy* spill)
  {
    return spill->slc;
  }

  const caf::Proxy<std::vector<caf::SRTrack>>&
  GetTracks(const caf::SRSliceRecoBranchProxy* reco)
  {
    return reco->trk;
  }

  const caf::Proxy<std::vector<caf::SRShower>>&
  GetShowers(const caf::SRSliceRecoBranchProxy* reco)
  {
    return reco->shw;
  }

  const caf::Proxy<std::vector<caf::SRStub>>&
  GetStubs(const caf::SRSliceRecoBranchProxy* reco)
  {
    return reco->stub;
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  // Instantiations
  template class VectorAdaptor<caf::StandardRecord, caf::SRSlice>;

  template class VectorAdaptor<caf::SRSliceRecoBranch, caf::SRTrack>;
  template class VectorAdaptor<caf::SRSliceRecoBranch, caf::SRShower>;
  template class VectorAdaptor<caf::SRSliceRecoBranch, caf::SRStub>;

  template class EnsembleVectorAdaptor<caf::StandardRecord, caf::SRSlice>;

  template class EnsembleVectorAdaptor<caf::SRSliceRecoBranch, caf::SRTrack>;
  template class EnsembleVectorAdaptor<caf::SRSliceRecoBranch, caf::SRShower>;
  template class EnsembleVectorAdaptor<caf::SRSliceRecoBranch, caf::SRStub>;
}
