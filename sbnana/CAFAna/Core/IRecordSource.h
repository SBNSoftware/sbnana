#pragma once

#include "CAFAnaCore/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/IRecordSink.h"

#include "sbnana/CAFAna/Core/SystShifts.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  using ISpillSource = beta::_IRecordSource<caf::SRSpillProxy>;
  using ISliceSource = beta::_IRecordSource<caf::SRSliceProxy>;

  using ISpillEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSpillProxy>;
  using ISliceEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSliceProxy>;

  //----------------------------------------------------------------------

  template<> class beta::_IRecordSource<caf::SRSliceProxy> : public beta::_IRecordSourceDefaultImpl<caf::SRSliceProxy>
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

  class SliceAdaptor : public beta::PassthroughExposure<ISpillSink, ISliceSource>
  {
  public:
    SliceAdaptor(ISpillSource& src);

    virtual void HandleRecord(const caf::SRSpillProxy* spill,
                              double weight) override;

    // Will need an EnsembleSliceAdaptor if we ever have syst weights that
    // apply at the spill level.
  };


  // Spill sources also provide a slice source (which loops over the slices)
  template<> class beta::_IRecordSource<caf::SRSpillProxy> : public beta::_IRecordSourceDefaultImpl<caf::SRSpillProxy>
  {
  public:
    ISliceSource& Slices() {return fSlices;}

  protected:
    _IRecordSource() : fSlices(*this) {}

    SliceAdaptor fSlices;
  };

  //----------------------------------------------------------------------

  // Introduce some aliases so we can express ourselves more succinctly

  using ITrackSource = beta::_IRecordSource<caf::SRTrackProxy>;
  using IShowerSource = beta::_IRecordSource<caf::SRShowerProxy>;
  using IStubSource = beta::_IRecordSource<caf::SRStubProxy>;

  using ISliceRecoBranchSource = beta::_IRecordSource<caf::SRSliceRecoBranchProxy>;
  using ISliceRecoBranchSink = beta::_IRecordSink<caf::SRSliceRecoBranchProxy>;

  //----------------------------------------------------------------------

  // Provide ability to get track / shower / stub sources from the reco branch
  // source.

  class TrackAdaptor : public beta::PassthroughExposure<ISliceRecoBranchSink,
                                                        ITrackSource>
  {
  public:
    TrackAdaptor(ISliceRecoBranchSource& src);

    virtual void HandleRecord(const caf::SRSliceRecoBranchProxy* reco,
                              double weight) override;
  };

  class ShowerAdaptor : public beta::PassthroughExposure<ISliceRecoBranchSink,
                                                         IShowerSource>
  {
  public:
    ShowerAdaptor(ISliceRecoBranchSource& src);

    virtual void HandleRecord(const caf::SRSliceRecoBranchProxy* reco,
                              double weight) override;
  };

  class StubAdaptor : public beta::PassthroughExposure<ISliceRecoBranchSink,
                                                       IStubSource>
  {
  public:
    StubAdaptor(ISliceRecoBranchSource& src);

    virtual void HandleRecord(const caf::SRSliceRecoBranchProxy* reco,
                              double weight) override;
  };


  template<> class beta::_IRecordSource<caf::SRSliceRecoBranchProxy> : public beta::_IRecordSourceDefaultImpl<caf::SRSliceRecoBranchProxy>
  {
  public:
    ITrackSource& Tracks() {return fTracks;}
    IShowerSource& Showers() {return fShowers;}
    IStubSource& Stubs() {return fStubs;}

  protected:
    _IRecordSource() : fTracks(*this), fShowers(*this), fStubs(*this) {}

    TrackAdaptor fTracks;
    ShowerAdaptor fShowers;
    StubAdaptor fStubs;
  };

  //----------------------------------------------------------------------

  using ITrackEnsembleSource = beta::_IRecordEnsembleSource<caf::SRTrackProxy>;
  using IShowerEnsembleSource = beta::_IRecordEnsembleSource<caf::SRShowerProxy>;
  using IStubEnsembleSource = beta::_IRecordEnsembleSource<caf::SRStubProxy>;

  using ISliceRecoBranchEnsembleSource = beta::_IRecordEnsembleSource<caf::SRSliceRecoBranchProxy>;
  using ISliceRecoBranchEnsembleSink = beta::_IRecordEnsembleSink<caf::SRSliceRecoBranchProxy>;

  //----------------------------------------------------------------------

  class TrackEnsembleAdaptor
    : public beta::PassthroughExposure<ISliceRecoBranchEnsembleSink,
                                       ITrackEnsembleSource>
  {
  public:
    TrackEnsembleAdaptor(ISliceRecoBranchEnsembleSource& src);

    virtual void HandleSingleRecord(const caf::SRSliceRecoBranchProxy* reco,
                                    double weight,
                                    int universeIdx) override;

    virtual void HandleEnsemble(const caf::SRSliceRecoBranchProxy* reco,
                                const std::vector<double>& weights) override;
  };

  class ShowerEnsembleAdaptor
    : public beta::PassthroughExposure<ISliceRecoBranchEnsembleSink,
                                       IShowerEnsembleSource>
  {
  public:
    ShowerEnsembleAdaptor(ISliceRecoBranchEnsembleSource& src);

    virtual void HandleSingleRecord(const caf::SRSliceRecoBranchProxy* reco,
                                    double weight,
                                    int universeIdx) override;

    virtual void HandleEnsemble(const caf::SRSliceRecoBranchProxy* reco,
                                const std::vector<double>& weights) override;
  };

  class StubEnsembleAdaptor
    : public beta::PassthroughExposure<ISliceRecoBranchEnsembleSink,
                                       IStubEnsembleSource>
  {
  public:
    StubEnsembleAdaptor(ISliceRecoBranchEnsembleSource& src);

    virtual void HandleSingleRecord(const caf::SRSliceRecoBranchProxy* reco,
                                    double weight,
                                    int universeIdx) override;

    virtual void HandleEnsemble(const caf::SRSliceRecoBranchProxy* reco,
                                const std::vector<double>& weights) override;
  };


  template<> class beta::_IRecordEnsembleSource<caf::SRSliceRecoBranchProxy> : public beta::_IRecordEnsembleSourceDefaultImpl<caf::SRSliceRecoBranchProxy>
  {
  public:
    ITrackEnsembleSource& Tracks() {return fTracks;}
    IShowerEnsembleSource& Showers() {return fShowers;}
    IStubEnsembleSource& Stubs() {return fStubs;}

  protected:
    _IRecordEnsembleSource() : fTracks(*this), fShowers(*this), fStubs(*this) {}

    TrackEnsembleAdaptor fTracks;
    ShowerEnsembleAdaptor fShowers;
    StubEnsembleAdaptor fStubs;
  };
}
