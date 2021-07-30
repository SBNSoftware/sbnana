#include "sbnana/CAFAna/Core/Loaders.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include <cassert>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  template<class SrcT> Sources<SrcT>::Sources()
  {
  }

  //----------------------------------------------------------------------
  template<class SrcT> Sources<SrcT>::~Sources()
  {
    // for(auto it: fSources) delete it.second;
  }

  //----------------------------------------------------------------------
  template<class SrcT> void Sources<SrcT>::SetLoaderPath(const std::string& path,
                                                         DataMC datamc,
                                                         SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Clear out the old one if necessary
    DisableLoader(datamc, swap);

    fLoaderPaths[key] = path;
  }

  //----------------------------------------------------------------------
  template<class SrcT> void Sources<SrcT>::SetLoaderFiles(const std::vector<std::string>& files,
                                                          DataMC datamc,
                                                          SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Clear out the old one if necessary
    DisableLoader(datamc, swap);

    fLoaderFiles[key] = files;
  }

  //----------------------------------------------------------------------
  template<class SrcT> void Sources<SrcT>::AddLoader(SrcT* src,
                                                     DataMC datamc,
                                                     SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Clear out the old one if necessary
    DisableLoader(datamc, swap);

    fSources[key] = src;
  }

  //----------------------------------------------------------------------
  template<class SrcT> void Sources<SrcT>::DisableLoader(DataMC datamc,
                                                         SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Clear out the current one if possible
    auto it = fSources.find(key);
    if(it != fSources.end()){
      delete it->second;
      fSources.erase(it);
    }

    fLoaderPaths.erase(key);
    fLoaderFiles.erase(key);
  }

  //----------------------------------------------------------------------
  template<class SrcT> SrcT& Sources<SrcT>::GetLoader(DataMC datamc,
                                                      SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Look up and return. Use kNullLoader if no loader is set for this config
    auto itLoader = fSources.find(key);
    if(itLoader != fSources.end()) return *itLoader->second;

    auto itPath = fLoaderPaths.find(key);
    if(itPath != fLoaderPaths.end()){
      fSources[key] = new SpectrumLoader(itPath->second);
      return *fSources[key];
    }
    auto itFiles = fLoaderFiles.find(key);
    if(itFiles != fLoaderFiles.end()){
      fSources[key] = new SpectrumLoader(itFiles->second);
      return *fSources[key];
    }

    return kNullLoader;
  }

  // Instantiate the ones we need
  template class Sources<SpectrumLoader>;
  template class Sources<SBNSpillSource>;
  template class Sources<ISpillSource>;
  template class Sources<ISliceSource>;
}
