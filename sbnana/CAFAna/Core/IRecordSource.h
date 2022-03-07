#pragma once

#include "CAFAnaCore/CAFAna/Core/IRecordSource.h"

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

  // Spill sources also provide a slice source (which loops over the slices)
  template<> class beta::_IRecordSource<caf::SRSpillProxy>
    : public beta::_IRecordSourceDefaultImpl<caf::SRSpillProxy>
  {
  public:
    ISliceSource& Slices() {return *fSlices;}

  protected:
    _IRecordSource();

    std::unique_ptr<ISliceSource> fSlices;
  };

  //----------------------------------------------------------------------

  // Provide ability to get track / shower / stub sources from the reco branch
  // source.

  template<> class beta::_IRecordSource<caf::SRSliceRecoBranchProxy>
    : public beta::_IRecordSourceDefaultImpl<caf::SRSliceRecoBranchProxy>
  {
  public:
    ITrackSource& Tracks() {return *fTracks;}
    IShowerSource& Showers() {return *fShowers;}
    IStubSource& Stubs() {return *fStubs;}

  protected:
    _IRecordSource();

    std::unique_ptr<ITrackSource> fTracks;
    std::unique_ptr<IShowerSource> fShowers;
    std::unique_ptr<IStubSource> fStubs;
  };

  //----------------------------------------------------------------------

  // Provide ability to get track / shower / stub sources from the reco branch
  // ensemble source.

  template<> class beta::_IRecordEnsembleSource<caf::SRSliceRecoBranchProxy>
    : public beta::_IRecordEnsembleSourceDefaultImpl<caf::SRSliceRecoBranchProxy>
  {
  public:
    ITrackEnsembleSource& Tracks() {return *fTracks;}
    IShowerEnsembleSource& Showers() {return *fShowers;}
    IStubEnsembleSource& Stubs() {return *fStubs;}

  protected:
    _IRecordEnsembleSource();

    std::unique_ptr<ITrackEnsembleSource> fTracks;
    std::unique_ptr<IShowerEnsembleSource> fShowers;
    std::unique_ptr<IStubEnsembleSource> fStubs;
  };
}
