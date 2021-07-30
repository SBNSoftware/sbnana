#pragma once

#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include <map>

namespace ana
{
  class SpectrumLoader;

  // TODO should these be in some sort of namespace?
  enum DataMC{kData, kMC, kNumDataMCs};
  enum SwappingConfig{kNonSwap, kNueSwap, kNuTauSwap, kIntrinsic, kNumSwappingConfigs};

  /// \brief Collection of SpectrumLoaders for many configurations
  template<class SrcT> class Sources
  {
  public:
    /// No loaders initialized. Use \ref SetLoaderPath to configure
    Sources();
    ~Sources();

    /// Configure loader via wildcard \a path
    void SetLoaderPath(const std::string& path,
                       DataMC datamc,
                       SwappingConfig swap = kNonSwap);

    /// Configure loader via explicit file list
    void SetLoaderFiles(const std::vector<std::string>& files,
                        DataMC datamc,
                        SwappingConfig swap = kNonSwap);

    void AddLoader(SrcT*,
                   DataMC datamc,
                   SwappingConfig swap = kNonSwap);

    void DisableLoader(DataMC datamc,
                       SwappingConfig swap = kNonSwap);

    /// Retrieve a specific loader
    SrcT& GetLoader(DataMC datamc,
                    SwappingConfig swap = kNonSwap);

    template<class T> auto& operator[](const T& x)
    {
      auto ret = new Sources<std::remove_reference_t<decltype((*fSources.begin()->second)[x])>>;

      for(int dmc = 0; dmc < kNumDataMCs; ++dmc){
        for(int swap = 0; swap < kNumSwappingConfigs; ++swap){
          if(dmc == kData && swap != kNonSwap) continue;
          ret->AddLoader(&GetLoader(DataMC(dmc), SwappingConfig(swap))[x],
                         DataMC(dmc), SwappingConfig(swap));
        }
      }

      return *ret;
    }

  protected:
    typedef std::tuple<DataMC, SwappingConfig> Key_t;

    // Hold a list of paths that have been set
    std::map<Key_t, std::string> fLoaderPaths;
    std::map<Key_t, std::vector<std::string>> fLoaderFiles;
    // Only reify them when someone actually calls GetLoader()
    std::map<Key_t, SrcT*> fSources;
  };


  using SpillSources = Sources<ISpillSource>;
  using SliceSources = Sources<ISliceSource>;


  class Loaders: public Sources<SpectrumLoader>
  {
  public:
    operator SliceSources&()
    {
      SliceSources* ret = new SliceSources;
      for(int dmc = 0; dmc < kNumDataMCs; ++dmc){
        for(int swap = 0; swap < kNumSwappingConfigs; ++swap){
          if(dmc == kData && swap != kNonSwap) continue;
          ret->AddLoader(&GetLoader(DataMC(dmc), SwappingConfig(swap)),
                         DataMC(dmc), SwappingConfig(swap));
        }
      }
      return *ret;
    }

    /// Call Go() on all the loaders
    void Go()
    {
      for(auto it: fSources) it.second->Go();
    }
  };

} // namespace
