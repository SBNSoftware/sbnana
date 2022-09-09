#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"

#include "cafanacore/Progress.h"
#include "cafanacore/ReweightableSpectrum.h"
#include "cafanacore/SAMQuerySource.h"
#include "cafanacore/SAMProjectSource.h"
#include "cafanacore/Spectrum.h"
#include "cafanacore/WildcardSource.h"

#include "sbnana/CAFAna/Core/Utilities.h"

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
  SpectrumLoaderBase::SpectrumLoaderBase()
    : fGone(false), fPOT(0), fPOTFromHist(0), fNReadouts(0)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::string& wildcard)
    : SpectrumLoaderBase()
  {
    fWildcard = wildcard;
    fFileSource = std::unique_ptr<IFileSource>(WildcardOrSAMQuery(wildcard));
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::SpectrumLoaderBase(const std::vector<std::string>& fnames) : SpectrumLoaderBase()
  {
    fWildcard = "file list";
    fFileSource = std::unique_ptr<IFileSource>(new FileListSource(fnames));

    assert(!fnames.empty());
    std::cout << "Loading from " << fnames.size() << " files" << std::endl;
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase::~SpectrumLoaderBase()
  {
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
} // namespace
