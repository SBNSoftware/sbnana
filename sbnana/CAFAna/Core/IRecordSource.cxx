#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/Multiverse.h"

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
                               const Multiverse& multiverse);

    virtual void HandleRecord(const caf::SRSliceProxy* slc, double weight) override;

    virtual const ana::FitMultiverse* GetMultiverse() const override {return fMultiverse;}

  protected:
    const Multiverse* fMultiverse;
    std::vector<SystShifts> fShifts;
  };

  //----------------------------------------------------------------------
  ShiftedSliceEnsembleSource::
  ShiftedSliceEnsembleSource(ISliceSource& src, const Multiverse& multiverse)
    : fMultiverse(&multiverse)
  {
    src.Register(this);

    // Turn the universes into concrete SystShifts objects up-front
    fShifts.reserve(multiverse.NUniv());
    for(unsigned int i = 0; i < multiverse.NUniv(); ++i) fShifts.emplace_back(multiverse.GetUniverse(i));
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
  Ensemble(const Multiverse& multiverse)
  {
    return fEnsembleSources.template Get<ShiftedSliceEnsembleSource>(&multiverse, *this, multiverse);
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
    : fSource(&src), fVecGetter(vecGetter)
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

  const caf::Proxy<std::vector<caf::SRTrueInteraction>>&
  GetNuTruths(const caf::SRSpillProxy* spill)
  {
    return spill->mc.nu;
  }

  const caf::Proxy<std::vector<caf::SRTrack>>&
  GetTracks(const caf::SRSliceProxy* slc)
  {
    return slc->reco.trk;
  }

  const caf::Proxy<std::vector<caf::SRShower>>&
  GetShowers(const caf::SRSliceProxy* slc)
  {
    return slc->reco.shw;
  }

  const caf::Proxy<std::vector<caf::SRStub>>&
  GetStubs(const caf::SRSliceProxy* slc)
  {
    return slc->reco.stub;
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  // Instantiations
  template class VectorAdaptor<caf::StandardRecord, caf::SRSlice>;
  template class VectorAdaptor<caf::StandardRecord, caf::SRTrueInteraction>;

  template class VectorAdaptor<caf::SRSlice, caf::SRTrack>;
  template class VectorAdaptor<caf::SRSlice, caf::SRShower>;
  template class VectorAdaptor<caf::SRSlice, caf::SRStub>;

  template class EnsembleVectorAdaptor<caf::StandardRecord, caf::SRSlice>;

  template class EnsembleVectorAdaptor<caf::SRSlice, caf::SRTrack>;
  template class EnsembleVectorAdaptor<caf::SRSlice, caf::SRShower>;
  template class EnsembleVectorAdaptor<caf::SRSlice, caf::SRStub>;
}
