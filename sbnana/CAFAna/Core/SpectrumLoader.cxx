#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/ReweightableSpectrum.h"
#include "CAFAna/Core/SAMProjectSource.h"
#include "CAFAna/Core/Spectrum.h"

#include "sbnana/CAFAna/Core/IRecordSink.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

namespace ana
{
  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::string& wildcard, int max)
    : SpectrumLoaderBase(wildcard), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::vector<std::string>& fnames, int max)
    : SpectrumLoaderBase(fnames), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader()
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader* SpectrumLoader::FromSAMProject(const std::string& proj,
                                                 int fileLimit)
  {
    SpectrumLoader* ret = new SpectrumLoader;
    ret->fWildcard = "project "+proj;
    ret->fFileSource = std::unique_ptr<IFileSource>(new SAMProjectSource(proj, fileLimit));
    return ret;
  }

  //----------------------------------------------------------------------
  SpectrumLoader::~SpectrumLoader()
  {
  }

  /*
  struct CompareByID
  {
    bool operator()(const Cut& a, const Cut& b) const
    {
      return a.ID() < b.ID();
    }
  };
  */

  //----------------------------------------------------------------------
  void SpectrumLoader::Go()
  {
    if(fGone){
      std::cerr << "Error: can only call Go() once on a SpectrumLoader" << std::endl;
      abort();
    }
    fGone = true;

    // Find all the unique cuts
    //    std::set<Cut, CompareByID> cuts;
    //    for(auto& shiftdef: fHistDefs)
    //      for(auto& cutdef: shiftdef.second)
    //        cuts.insert(cutdef.first);
    //    for(const Cut& cut: cuts) fAllCuts.push_back(cut);

    //    fLivetimeByCut.resize(fAllCuts.size());
    //    fPOTByCut.resize(fAllCuts.size());


    const int Nfiles = NFiles();

    Progress* prog = 0;

    caf::SRBranchRegistry::clear();

    int fileIdx = -1;
    while(TFile* f = GetNextFile()){
      ++fileIdx;

      if(Nfiles >= 0 && !prog){
        int totsinks = 0;
        for(const ISpillSink* s: ISpillSource::fSinks) totsinks += s->NSinks();
        for(const ISliceSink* s: ISliceSource::fSinks) totsinks += s->NSinks();

        prog = new Progress(TString::Format("Filling %d spectra from %d files matching '%s'", totsinks, Nfiles, fWildcard.c_str()).Data());
      }

      HandleFile(f, Nfiles == 1 ? prog : 0);

      if(Nfiles > 1 && prog) prog->SetProgress((fileIdx+1.)/Nfiles);
    } // end for fileIdx

    if(prog){
      prog->Done();
      delete prog;
    }

    StoreExposures(); // also triggers the POT printout

    // TODO TODO any need for a RemoveLoader() call now?
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleFile(TFile* f, Progress* prog)
  {
    assert(!f->IsZombie());

    // Test for multi (trees are inside a directory) or single (tree is at top
    // level) tree cases.
    TDirectory* dir = 0;
    TTree* tr = 0;

    TObject* obj = f->Get("recTree");
    assert(obj); // Must exist in one form or the other

    // It might seem like you could use GetObject() to do the type-checking
    // directly, but that method seems to have a memory leak in the case that
    // the types don't match.
    if(obj->ClassName() == std::string("TTree")){
      tr = (TTree*)obj;
    }
    else{
      dir = (TDirectory*)obj;
      tr = (TTree*)dir->Get("rec");
      assert(tr);
    }

    const caf::CAFType type = caf::GetCAFType(dir, tr);

    long n;
    caf::SRSpillProxy sr(dir, tr, "rec", n, 0);

    //    FloatingExceptionOnNaN fpnan;

    long Nentries = tr->GetEntries();
    if (max_entries != 0 && max_entries < Nentries) Nentries = max_entries;

    for(n = 0; n < Nentries; ++n){
      if(type != caf::kFlatMultiTree) tr->LoadTree(n); // for all single-tree modes

      HandleRecord(&sr);

      if(prog) prog->SetProgress(double(n)/Nentries);
    } // end for n
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleRecord(caf::SRSpillProxy* sr)
  {
    if(sr->hdr.first_in_subrun){
      fPOT += sr->hdr.pot;
      // TODO think about if this should be gated behind first_in_file. At the
      // moment I think these will be synonymous. And despite the comment on
      // hdr.pot, I think it may be file-based in practice too.
      fNGenEvt += sr->hdr.ngenevt;
    }

    for(ISpillSink* s: ISpillSource::fSinks) s->HandleRecord(sr, 1);

    for(caf::SRSliceProxy& slc: sr->slc){
      for(ISliceSink* s: ISliceSource::fSinks){
        s->HandleRecord(&slc, 1);
      }
    }
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::StoreExposures()
  {
    if(fabs(fPOT - fPOTFromHist)/std::min(fPOT, fPOTFromHist) > 0.001){
      std::cout << fPOT << " POT from hdr differs from " << fPOTFromHist << " POT from the TotalPOT histogram!" << std::endl;
      abort();
    }

    std::cout << fPOT << " POT over " << fNGenEvt << " readouts" << std::endl;


    for(ISpillSink* s: ISpillSource::fSinks){
      s->HandlePOT(fPOT);
      s->HandleLivetime(fNGenEvt);
    }

    for(ISliceSink* s: ISliceSource::fSinks){
      s->HandlePOT(fPOT);
      s->HandleLivetime(fNGenEvt);
    }
  }

  //----------------------------------------------------------------------
  NullLoader::NullLoader()
  {
  }

  //----------------------------------------------------------------------
  NullLoader::~NullLoader()
  {
  }

  //----------------------------------------------------------------------
  void NullLoader::Go()
  {
  }

} // namespace
