#include "sbnana/CAFAna/Core/Loaders.h"

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include <cassert>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  Loaders::Loaders()
  {
  }

  //----------------------------------------------------------------------
  Loaders::~Loaders()
  {
    // for(auto it: fLoaders) delete it.second;
  }

  //----------------------------------------------------------------------
  void Loaders::SetLoaderPath(const std::string& path,
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
  void Loaders::SetLoaderFiles(const std::vector<std::string>& files,
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
  void Loaders::AddLoader(SpectrumLoaderBase* file,
                               DataMC datamc,
                               SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Clear out the old one if necessary
    DisableLoader(datamc, swap);

    fLoaders[key] = file;
  }

  //----------------------------------------------------------------------
  void Loaders::DisableLoader(DataMC datamc,
                              SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Clear out the current one if possible
    auto it = fLoaders.find(key);
    if(it != fLoaders.end()){
      delete it->second;
      fLoaders.erase(it);
    }

    fLoaderPaths.erase(key);
    fLoaderFiles.erase(key);
  }

  //----------------------------------------------------------------------
  SpectrumLoaderBase& Loaders::GetLoader(DataMC datamc,
                                         SwappingConfig swap)
  {
    assert(datamc == kMC || swap == kNonSwap);

    const Key_t key(datamc, swap);

    // Look up and return. Use fNull if no loader is set for this config
    auto itLoader = fLoaders.find(key);
    if(itLoader != fLoaders.end()) return *itLoader->second;

    auto itPath = fLoaderPaths.find(key);
    if(itPath != fLoaderPaths.end()){
      fLoaders[key] = new SpectrumLoader(itPath->second);
      return *fLoaders[key];
    }
    auto itFiles = fLoaderFiles.find(key);
    if(itFiles != fLoaderFiles.end()){
      fLoaders[key] = new SpectrumLoader(itFiles->second);
      return *fLoaders[key];
    }

    return fNull;
  }

  //----------------------------------------------------------------------
  void Loaders::Go()
  {
    for(auto it: fLoaders) it.second->Go();
  }
}
