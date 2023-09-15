#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"

#include "sbnana/CAFAna/Core/Progress.h"
#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"
#include "sbnana/CAFAna/Core/SAMQuerySource.h"
#include "sbnana/CAFAna/Core/SAMProjectSource.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/WildcardSource.h"
#include "sbnana/CAFAna/Core/Tree.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "ifdh.h"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include "boost/algorithm/string.hpp"

namespace ana
{
  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::Erase(Spectrum* s)
  {
    auto it = std::find(spects.begin(), spects.end(), s);
    if(it != spects.end()) spects.erase(it);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::Erase(ReweightableSpectrum* rs)
  {
    auto it = std::find(rwSpects.begin(), rwSpects.end(), rs);
    if(it != rwSpects.end()) rwSpects.erase(it);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::RemoveLoader(SpectrumLoaderBase* l)
  {
    for(Spectrum* s: spects) s->RemoveLoader(l);
    for(ReweightableSpectrum* rs: rwSpects) rs->RemoveLoader(l);
  }

  //----------------------------------------------------------------------
  size_t SpectrumLoaderBase::SpectList::TotalSize() const
  {
    return spects.size() + rwSpects.size();
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::GetSpectra(std::vector<Spectrum*>& ss)
  {
    ss.insert(ss.end(), spects.begin(), spects.end());
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::SpectList::
  GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss)
  {
    ss.insert(ss.end(), rwSpects.begin(), rwSpects.end());
  }

  //----------------------------------------------------------------------
  template<class T, class U> U& SpectrumLoaderBase::IDMap<T, U>::
  operator[](const T& key)
  {
    for(auto& it: fElems){
      if(it.first.ID() == key.ID()) return it.second;
    }
    fElems.push_back(std::make_pair(key, U()));
    return fElems.back().second;
  }

  //----------------------------------------------------------------------
  template<class T, class U> template<class V> void SpectrumLoaderBase::IDMap<T, U>::Erase(const V& v)
  {
    for(auto& it: fElems) it.second.Erase(v);
  }

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::
  RemoveLoader(SpectrumLoaderBase* l)
  {
    for(auto& it: fElems) it.second.RemoveLoader(l);
  }

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::Clear()
  {
    fElems.clear();
  }

  //----------------------------------------------------------------------
  template<class T, class U> size_t SpectrumLoaderBase::IDMap<T, U>::
  TotalSize()
  {
    size_t ret = 0;
    for(auto& it: fElems) ret += it.second.TotalSize();
    return ret;
  }

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::
  GetSpectra(std::vector<Spectrum*>& ss)
  {
    for(auto& it: fElems) it.second.GetSpectra(ss);
  }

  //----------------------------------------------------------------------
  template<class T, class U> void SpectrumLoaderBase::IDMap<T, U>::
  GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss)
  {
    for(auto& it: fElems) it.second.GetReweightableSpectra(ss);
  }

  // Start of SpectrumLoaderBase proper

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(DataSource src)
    : fSource(src), fGone(false), fPOT(0), fPOTFromHist(0), fNReadouts(0)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::string& wildcard,
                                         DataSource src)
    : SpectrumLoaderBase(src)
  {
    fWildcard = wildcard;
    fFileSource = std::unique_ptr<IFileSource>(WildcardOrSAMQuery(wildcard));
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::vector<std::string>& fnames,
                                         DataSource src)
    : SpectrumLoaderBase(src)
  {
    fWildcard = "file list";
    fFileSource = std::unique_ptr<IFileSource>(new FileListSource(fnames));

    assert(!fnames.empty());
    std::cout << "Loading from " << fnames.size() << " files" << std::endl;
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::~SpectrumLoaderBase()
  {
    fHistDefs.RemoveLoader(this);
  }

  //----------------------------------------------------------------------
  IFileSource* SpectrumLoaderBase::
  WildcardOrSAMQuery(const std::string& str) const
  {
    int stride = -1;
    int offset = -1;
    if(getenv("CAFANA_STRIDE")){
      stride = atoi(getenv("CAFANA_STRIDE"));
      if(stride > 1 && getenv("CAFANA_OFFSET")){
        offset = atoi(getenv("CAFANA_OFFSET"));
      }
    }

    // stat() blows up on strings with spaces
    if(str.find(' ') == std::string::npos){
      WildcardSource* ret = new WildcardSource(str, stride, offset);
      if(ret->NFiles() > 0) return ret;
      delete ret;
    }

    // Maybe this the name of a SAM project?
    {
      IFDHSilent silent; // the usual case is for this to fail
      ifdh_ns::ifdh i;
      // findProject always gives back an address just by gluing bits together.
      // (a tad annoying, because it _does_ go and look for the project, and
      // even would print its 'didn't-find-this-project' error out to stderr if
      // not for the IFDHSilent, but without scraping its stderr there's no way
      // to know whether the project is there or not -- you still get the URL.)
      // however, the WebAPI call looking for the /status will return a 404 if
      // the project doesn't exist.  (suggested by Robert I. in
      // INC000000925362)
      try{
        ifdh_util_ns::WebAPI webapi(i.findProject(str, SAMExperiment()) +  "/status");
        return new SAMProjectSource(str);
      }
      catch(ifdh_util_ns::WebAPIException &e){
        ;
      }
    }

    // Maybe this is a SAM dataset or query?
    return new SAMQuerySource(str, stride, offset);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const Var& var,
                                       const SpillCut& spillcut,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fHistDefs[spillcut][shift][cut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const MultiVar& var,
                                       const SpillCut& spillcut,
                                       const Cut& cut,
                                       const SystShifts& shift,
                                       const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fHistDefs[spillcut][shift][cut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const SpillVar& var,
                                       const SpillCut& cut,
                                       const SpillVar& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fSpillHistDefs[cut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const SpillMultiVar& var,
                                       const SpillCut& cut,
                                       const SpillVar& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fSpillHistDefs[cut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const TruthVar& var,
                                       const TruthCut truthcut,
                                       const SpillCut& spillcut,
                                       const Cut& cut, // loop over reco slices and see if any matched to this truth and pass "cut"
                                       const TruthVar& wei)
  { 
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fTruthHistDefs[spillcut][cut][truthcut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
    
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddSpectrum(Spectrum& spect,
                                       const TruthMultiVar& var,
                                       const TruthCut truthcut,
                                       const SpillCut& spillcut,
                                       const Cut& cut, // loop over reco slices and see if any matched to this truth and pass "cut"
                                       const TruthVar& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fTruthHistDefs[spillcut][cut][truthcut][wei][var].spects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending

  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::RemoveSpectrum(Spectrum* spect)
  {
    fHistDefs.Erase(spect);
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddReweightableSpectrum(ReweightableSpectrum& spect,
                                                   const Var& var,
                                                   const Cut& cut,
                                                   const SystShifts& shift,
                                                   const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fHistDefs[kNoSpillCut][shift][cut][wei][var].rwSpects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddReweightableSpectrum(ReweightableSpectrum& spect,
                                                   const Var& var,
                                                   const SpillCut& spillcut,
                                                   const SliceCut& slicecut,
                                                   const SystShifts& shift,
                                                   const Var& wei)
  {
    if(fGone){
      std::cerr << "Error: can't add Spectra after the call to Go()" << std::endl;
      abort();
    }

    fHistDefs[spillcut][shift][slicecut][wei][var].rwSpects.push_back(&spect);

    spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddTree(Tree& tree,
                                   const std::vector<std::string>& labels,
                                   const std::vector<Var>& vars,
                                   const SpillCut& spillcut,
                                   const Cut& cut,
                                   const SystShifts& shift)
  {
    if(fGone){
      std::cerr << "Error: can't add Tree after the call to Go()" << std::endl;
      abort();
    }

    for ( unsigned int idx=0; idx<labels.size(); ++idx ) {
      fTreeDefs[spillcut][shift][cut][&tree][ vars.at(idx) ] = labels.at(idx);
    }

    // TODO do we need to add/remove loaders?
    //spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddTree(Tree& tree,
                                   const std::vector<std::string>& labels,
                                   const std::vector<MultiVar>& vars,
                                   const SpillCut& spillcut,
                                   const Cut& cut,
                                   const SystShifts& shift)
  {
    if(fGone){
      std::cerr << "Error: can't add Tree after the call to Go()" << std::endl;
      abort();
    }

    for ( unsigned int idx=0; idx<labels.size(); ++idx ) {
      fTreeDefs[spillcut][shift][cut][&tree][ vars.at(idx) ] = labels.at(idx);
    }

    // TODO do we need to add/remove loaders?
    //spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddTree(Tree& tree,
                                   const std::vector<std::string>& labels,
                                   const std::vector<SpillVar>& vars,
                                   const SpillCut& spillcut)
  {
    if(fGone){
      std::cerr << "Error: can't add Tree after the call to Go()" << std::endl;
      abort();
    }

    for ( unsigned int idx=0; idx<labels.size(); ++idx ) {
      fSpillTreeDefs[spillcut][&tree][ vars.at(idx) ] = labels.at(idx);
    }

    // TODO do we need to add/remove loaders?
    //spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddTree(Tree& tree,
                                   const std::vector<std::string>& labels,
                                   const std::vector<SpillMultiVar>& vars,
                                   const SpillCut& spillcut)
  {
    if(fGone){
      std::cerr << "Error: can't add Tree after the call to Go()" << std::endl;
      abort();
    }

    for ( unsigned int idx=0; idx<labels.size(); ++idx ) {
      fSpillTreeDefs[spillcut][&tree][ vars.at(idx) ] = labels.at(idx);
    }

    // TODO do we need to add/remove loaders?
    //spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddNSigmasTree(NSigmasTree& tree,
                                          const std::vector<std::string>& labels,
                                          const std::vector<const ISyst*>& systsToStore,
                                          const SpillCut& spillcut,
                                          const Cut& cut,
                                          const SystShifts& shift)
  {
    if(fGone){
      std::cerr << "Error: can't add NSigmasTree after the call to Go()" << std::endl;
      abort();
    }

    for ( unsigned int idx=0; idx<labels.size(); ++idx ) {
      fNSigmasTreeDefs[spillcut][shift][cut][&tree][ systsToStore.at(idx) ] = labels.at(idx);
    }

    // TODO do we need to add/remove loaders?
    //spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::AddNUniversesTree(NUniversesTree& tree,
                                             const std::vector<std::string>& labels,
                                             const std::vector<std::vector<Var>>& univKnobs,
                                             const SpillCut& spillcut,
                                             const Cut& cut,
                                             const SystShifts& shift)
  {
    if(fGone){
      std::cerr << "Error: can't add NUniversesTree after the call to Go()" << std::endl;
      abort();
    }

    for ( unsigned int idx=0; idx<labels.size(); ++idx ) {
      std::vector<VarOrMultiVar> vecVarOrMultiVar;
      for ( unsigned int idxUniv=0; idxUniv<univKnobs.at(idx).size(); ++idxUniv ) vecVarOrMultiVar.push_back( univKnobs.at(idx).at(idxUniv) );
      fNUniversesTreeDefs[spillcut][shift][cut][&tree][ vecVarOrMultiVar ] = labels.at(idx);
    }

    // TODO do we need to add/remove loaders?
    //spect.AddLoader(this); // Remember we have a Go() pending
  }

  //----------------------------------------------------------------------
  void SpectrumLoaderBase::
  RemoveReweightableSpectrum(ReweightableSpectrum* spect)
  {
    fHistDefs.Erase(spect);
  }

  //----------------------------------------------------------------------
  int SpectrumLoaderBase::NFiles() const
  {
    return fFileSource->NFiles();
  }

  //----------------------------------------------------------------------
  TFile* SpectrumLoaderBase::GetNextFile()
  {
    TFile* f = fFileSource->GetNextFile();
    if(!f) return 0; // out of files

    TH1* hPOT = (TH1*)f->Get("TotalPOT");
    assert(hPOT);
    fPOTFromHist  += hPOT->Integral(0, -1);

    return f;
  }

  //----------------------------------------------------------------------
  void NullLoader::Go()
  {
  }

  //----------------------------------------------------------------------
  NullLoader::~NullLoader()
  {
  }

  // Apparently the existence of fSpillDefs isn't enough and I need to spell
  // this out to make sure the function bodies are generated.
  template struct SpectrumLoaderBase::IDMap<SpillCut, SpectrumLoaderBase::IDMap<SystShifts, SpectrumLoaderBase::IDMap<Cut, SpectrumLoaderBase::IDMap<Var, SpectrumLoaderBase::IDMap<SpectrumLoaderBase::VarOrMultiVar, SpectrumLoaderBase::SpectList>>>>>;

  template struct SpectrumLoaderBase::IDMap<SpillCut, SpectrumLoaderBase::IDMap<SpillVar, SpectrumLoaderBase::IDMap<SpectrumLoaderBase::SpillVarOrMultiVar, SpectrumLoaderBase::SpectList>>>;
} // namespace
