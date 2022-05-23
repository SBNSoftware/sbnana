#pragma once

#include "cafanacore/IRecordSource.h"

#include "sbnana/CAFAna/Core/IRecordSink.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  class ISyst;
  template<class SystT> class _Multiverse;
  using Multiverse = _Multiverse<ISyst>;

  class SystShifts;

  //----------------------------------------------------------------------
  // Introduce some aliases so we can express ourselves more succinctly

  using ISpillSource = beta::_IRecordSource<caf::SRSpillProxy>;
  using ISliceSource = beta::_IRecordSource<caf::SRSliceProxy>;
  using INuTruthSource = beta::_IRecordSource<caf::SRTrueInteractionProxy>;

  using ISpillEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSpillProxy>;
  using ISliceEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSliceProxy>;
  using INuTruthEnsembleSource = beta::_IRecordEnsembleSource<caf::SRTrueInteractionProxy>;

  //----------------------------------------------------------------------

  using ITrackSource = beta::_IRecordSource<caf::SRTrackProxy>;
  using IShowerSource = beta::_IRecordSource<caf::SRShowerProxy>;
  using IStubSource = beta::_IRecordSource<caf::SRStubProxy>;

  //----------------------------------------------------------------------

  using ITrackEnsembleSource = beta::_IRecordEnsembleSource<caf::SRTrackProxy>;
  using IShowerEnsembleSource = beta::_IRecordEnsembleSource<caf::SRShowerProxy>;
  using IStubEnsembleSource = beta::_IRecordEnsembleSource<caf::SRStubProxy>;

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

  const caf::Proxy<std::vector<caf::SRTrueInteraction>>& GetNuTruths(const caf::SRSpillProxy* spill);

  const caf::Proxy<std::vector<caf::SRTrack>>& GetTracks(const caf::SRSliceProxy* reco);

  const caf::Proxy<std::vector<caf::SRShower>>& GetShowers(const caf::SRSliceProxy* reco);

  const caf::Proxy<std::vector<caf::SRStub>>& GetStubs(const caf::SRSliceProxy* reco);

  //----------------------------------------------------------------------

  template<> class beta::_IRecordSource<caf::SRSliceProxy>
    : public beta::_IRecordSourceDefaultImpl<caf::SRSliceProxy>
  {
  public:
    // Weight-based ensembles are still supported
    using _IRecordSourceDefaultImpl::Ensemble;

    // But also support an ensemble based on SystShifts
    ISliceEnsembleSource& Ensemble(const Multiverse& multiverse);

    ITrackSource& Tracks() {return fTracks;}
    IShowerSource& Showers() {return fShowers;}
    IStubSource& Stubs() {return fStubs;}

  protected:
    IDDict<const FitMultiverse*, ISliceEnsembleSource> fEnsembleSources;

    VectorAdaptor<caf::SRSlice, caf::SRTrack> fTracks{*this, GetTracks};
    VectorAdaptor<caf::SRSlice, caf::SRShower> fShowers{*this, GetShowers};
    VectorAdaptor<caf::SRSlice, caf::SRStub> fStubs{*this, GetStubs};
  };

  //----------------------------------------------------------------------

  // Spill sources also provide a slice source (which loops over the slices)
  template<> class beta::_IRecordSource<caf::SRSpillProxy>
    : public beta::_IRecordSourceDefaultImpl<caf::SRSpillProxy>
  {
  public:
    ISliceSource& Slices() {return fSlices;}
    INuTruthSource& NuTruths() {return fNuTruths;}

  protected:
    VectorAdaptor<caf::StandardRecord, caf::SRSlice> fSlices{*this, GetSlices};
    VectorAdaptor<caf::StandardRecord, caf::SRTrueInteraction> fNuTruths{*this, GetNuTruths};
  };
  //----------------------------------------------------------------------


  /// Helper class for implementing looping over slices, tracks, etc
  template<class FromT, class ToT> class EnsembleVectorAdaptor:
    public beta::PassthroughExposure<beta::_IRecordEnsembleSink<caf::Proxy<FromT>>,
                                     beta::_IRecordEnsembleSource<caf::Proxy<ToT>>>
  {
  public:
    using Source_t = beta::_IRecordEnsembleSource<caf::Proxy<FromT>>;
    using Func_t = std::function<const caf::Proxy<std::vector<ToT>>&(const caf::Proxy<FromT>*)>;

    EnsembleVectorAdaptor(Source_t& src, Func_t vecGetter);

    virtual void HandleSingleRecord(const caf::Proxy<FromT>* rec,
                                    double weight,
                                    int universeIdx) override;

    virtual void HandleEnsemble(const caf::Proxy<FromT>* rec,
                                const std::vector<double>& weight) override;

    virtual const ana::FitMultiverse* GetMultiverse() const {return fSource->GetMultiverse();}

  protected:
    const Source_t* fSource;
    Func_t fVecGetter;
  };

  //----------------------------------------------------------------------

  // Provide ability to get track / shower / stub sources from the reco branch
  // ensemble source.

  template<> class beta::_IRecordEnsembleSource<caf::SRSliceProxy>
    : public beta::_IRecordEnsembleSourceDefaultImpl<caf::SRSliceProxy>
  {
  public:
    ITrackEnsembleSource& Tracks() {return fTracks;}
    IShowerEnsembleSource& Showers() {return fShowers;}
    IStubEnsembleSource& Stubs() {return fStubs;}

  protected:
    EnsembleVectorAdaptor<caf::SRSlice, caf::SRTrack> fTracks{*this, GetTracks};
    EnsembleVectorAdaptor<caf::SRSlice, caf::SRShower> fShowers{*this, GetShowers};
    EnsembleVectorAdaptor<caf::SRSlice, caf::SRStub> fStubs{*this, GetStubs};
  };

  //----------------------------------------------------------------------
  // Spill sources also provide a slice source (which loops over the slices)
  template<> class beta::_IRecordEnsembleSource<caf::SRSpillProxy>
    : public beta::_IRecordEnsembleSourceDefaultImpl<caf::SRSpillProxy>
  {
  public:
    ISliceEnsembleSource& Slices() {return fSlices;}
    INuTruthEnsembleSource& NuTruths() {return fNuTruths;}

  protected:
    EnsembleVectorAdaptor<caf::StandardRecord, caf::SRSlice> fSlices{*this, GetSlices};
    EnsembleVectorAdaptor<caf::StandardRecord, caf::SRTrueInteraction> fNuTruths{*this, GetNuTruths};
  };
}
