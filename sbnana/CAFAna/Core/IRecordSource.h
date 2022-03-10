#pragma once

#include "cafanacore/IRecordSource.h"

#include "sbnana/CAFAna/Core/IRecordSink.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  class SystShifts;

  //----------------------------------------------------------------------
  // Introduce some aliases so we can express ourselves more succinctly

  using ISpillSource = beta::_IRecordSource<caf::SRSpillProxy>;
  using ISliceSource = beta::_IRecordSource<caf::SRSliceProxy>;

  using ISpillEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSpillProxy>;
  using ISliceEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSliceProxy>;

  //----------------------------------------------------------------------

  using ITrackSource = beta::_IRecordSource<caf::SRTrackProxy>;
  using IShowerSource = beta::_IRecordSource<caf::SRShowerProxy>;
  using IStubSource = beta::_IRecordSource<caf::SRStubProxy>;

  using ISliceRecoBranchSource = beta::_IRecordSource<caf::SRSliceRecoBranchProxy>;
  using ISliceRecoBranchSink = beta::_IRecordSink<caf::SRSliceRecoBranchProxy>;

  //----------------------------------------------------------------------

  using ITrackEnsembleSource = beta::_IRecordEnsembleSource<caf::SRTrackProxy>;
  using IShowerEnsembleSource = beta::_IRecordEnsembleSource<caf::SRShowerProxy>;
  using IStubEnsembleSource = beta::_IRecordEnsembleSource<caf::SRStubProxy>;

  using ISliceRecoBranchEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSliceRecoBranchProxy>;
  using ISliceRecoBranchEnsembleSink = beta::_IRecordEnsembleSink<caf::SRSliceRecoBranchProxy>;

  //----------------------------------------------------------------------

  template<> class beta::_IRecordSource<caf::SRSliceProxy>
    : public beta::_IRecordSourceDefaultImpl<caf::SRSliceProxy>
  {
  public:
    // Weight-based ensembles are still supported
    using _IRecordSourceDefaultImpl::Ensemble;

    // But also support an ensemble basd on SystShifts
    ISliceEnsembleSource& Ensemble(const std::vector<SystShifts>& shifts,
                                   int multiverseId);

  protected:
    IDDict<ISliceEnsembleSource> fEnsembleSources;
  };

  //----------------------------------------------------------------------

  /// Helper class for implementing looping over slices, tracks, etc
  template<class FromT, class ToT> class VectorAdaptor:
    public beta::PassthroughExposure<beta::_IRecordSink<caf::Proxy<FromT>>,
                                     beta::_IRecordSource<caf::Proxy<ToT>>>
  {
  public:
    typedef std::function<const caf::Proxy<std::vector<ToT>>&(const caf::Proxy<FromT>*)> Func_t;

    VectorAdaptor(beta::_IRecordSource<caf::Proxy<FromT>>& src, Func_t vecGetter);
    virtual void HandleRecord(const caf::Proxy<FromT>* rec, double weight) override;
  protected:
    Func_t fVecGetter;
  };

  //----------------------------------------------------------------------

  // Accessors needed by VectorAdaptor

  const caf::Proxy<std::vector<caf::SRSlice>>& GetSlices(const caf::SRSpillProxy* spill);

  const caf::Proxy<std::vector<caf::SRTrack>>& GetTracks(const caf::SRSliceRecoBranchProxy* reco);

  const caf::Proxy<std::vector<caf::SRShower>>& GetShowers(const caf::SRSliceRecoBranchProxy* reco);

  const caf::Proxy<std::vector<caf::SRStub>>& GetStubs(const caf::SRSliceRecoBranchProxy* reco);

  //----------------------------------------------------------------------

  // Spill sources also provide a slice source (which loops over the slices)
  template<> class beta::_IRecordSource<caf::SRSpillProxy>
    : public beta::_IRecordSourceDefaultImpl<caf::SRSpillProxy>
  {
  public:
    ISliceSource& Slices() {return fSlices;}

  protected:
    VectorAdaptor<caf::StandardRecord, caf::SRSlice> fSlices{*this, GetSlices};
  };

  //----------------------------------------------------------------------

  // Provide ability to get track / shower / stub sources from the reco branch
  // source.

  template<> class beta::_IRecordSource<caf::SRSliceRecoBranchProxy>
    : public beta::_IRecordSourceDefaultImpl<caf::SRSliceRecoBranchProxy>
  {
  public:
    ITrackSource& Tracks() {return fTracks;}
    IShowerSource& Showers() {return fShowers;}
    IStubSource& Stubs() {return fStubs;}

  protected:
    VectorAdaptor<caf::SRSliceRecoBranch, caf::SRTrack> fTracks{*this, GetTracks};
    VectorAdaptor<caf::SRSliceRecoBranch, caf::SRShower> fShowers{*this, GetShowers};
    VectorAdaptor<caf::SRSliceRecoBranch, caf::SRStub> fStubs{*this, GetStubs};
  };

  //----------------------------------------------------------------------

  /// Helper class for implementing looping over slices, tracks, etc
  template<class FromT, class ToT> class EnsembleVectorAdaptor:
    public beta::PassthroughExposure<beta::_IRecordEnsembleSink<caf::Proxy<FromT>>,
                                     beta::_IRecordEnsembleSource<caf::Proxy<ToT>>>
  {
  public:
    typedef std::function<const caf::Proxy<std::vector<ToT>>&(const caf::Proxy<FromT>*)> Func_t;

    EnsembleVectorAdaptor(beta::_IRecordEnsembleSource<caf::Proxy<FromT>>& src, Func_t vecGetter);

    virtual void HandleSingleRecord(const caf::Proxy<FromT>* rec,
                                    double weight,
                                    int universeIdx) override;

    virtual void HandleEnsemble(const caf::Proxy<FromT>* rec,
                                const std::vector<double>& weight) override;
  protected:
    Func_t fVecGetter;
  };

  //----------------------------------------------------------------------
  // Spill sources also provide a slice source (which loops over the slices)
  template<> class beta::_IRecordEnsembleSource<caf::SRSpillProxy>
    : public beta::_IRecordEnsembleSourceDefaultImpl<caf::SRSpillProxy>
  {
  public:
    ISliceEnsembleSource& Slices() {return fSlices;}

  protected:
    EnsembleVectorAdaptor<caf::StandardRecord, caf::SRSlice> fSlices{*this, GetSlices};
  };

  //----------------------------------------------------------------------

  // Provide ability to get track / shower / stub sources from the reco branch
  // ensemble source.

  template<> class beta::_IRecordEnsembleSource<caf::SRSliceRecoBranchProxy>
    : public beta::_IRecordEnsembleSourceDefaultImpl<caf::SRSliceRecoBranchProxy>
  {
  public:
    ITrackEnsembleSource& Tracks() {return fTracks;}
    IShowerEnsembleSource& Showers() {return fShowers;}
    IStubEnsembleSource& Stubs() {return fStubs;}

  protected:
    EnsembleVectorAdaptor<caf::SRSliceRecoBranch, caf::SRTrack> fTracks{*this, GetTracks};
    EnsembleVectorAdaptor<caf::SRSliceRecoBranch, caf::SRShower> fShowers{*this, GetShowers};
    EnsembleVectorAdaptor<caf::SRSliceRecoBranch, caf::SRStub> fStubs{*this, GetStubs};
  };
}
