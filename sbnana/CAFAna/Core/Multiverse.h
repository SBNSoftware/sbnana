#pragma once

#include "CAFAna/Core/INamed.h"
#include "CAFAna/Core/IFitSyst.h"

#include <unordered_map>
#include <vector>

class TDirectory;

namespace ana
{
  template<class SystT> using _Universe = std::unordered_map<const SystT*, double>;

  using FitUniverse = _Universe<IFitSyst>;

  /// \brief Collection of "universes" (SystShifts)
  ///
  /// Two Multiverses are equivalent if-and-only-if they have the same address
  ///
  /// Universe zero is always the nominal universe (no systematic shifts)
  class FitMultiverse: public INamed
  {
  public:
    // User should not delete these though
    virtual ~FitMultiverse() {}

    /// \brief Named constructor. Scan each parameter sequentially
    ///
    /// \param systs  The list of systematic parameters to scan
    /// \param nsigma Each parameter will be scanned from -nsigma to +nsigma
    static const FitMultiverse& Hypercross(const std::vector<const IFitSyst*>& systs,
                                           int nsigma = 3);

    static const unsigned int kTrulyRandom = 123456789;
    
    /// \brief Named constructor. Throw all parameters as gaussians
    ///
    /// \param systs The list of systematic parameters to vary
    /// \param Nuniv Number of universes to generate
    /// \param seed  Pseudo-random seed. Pass Multiverse::kTrulyRandom for a
    ///              real random seed. Beware incompatibility between universes
    ///              with different seeds.
    static const FitMultiverse& RandomGas(const std::vector<const IFitSyst*>& systs,
                                          int Nuniv,
                                          unsigned int seed);

    /// Total number of universes, including nominal at index 0
    unsigned int NUniv() const {return fUnivs.size();}

    /// Details of a particular universe
    const FitUniverse& GetUniverse(int i) const {return fUnivs[i];}

    void SaveTo(TDirectory* dir, const std::string& name) const;

    /// Usually these return unique_ptr, but Multiverses are globally managed
    static const FitMultiverse* LoadFrom(TDirectory* dir, const std::string& name);

  protected:
    FitMultiverse(const std::string& shortName,
                  const std::string& latexName,
                  const std::vector<FitUniverse>& univs);

    FitMultiverse(const FitMultiverse&) = delete;

    /// Special move constructor ONLY to help derived _Multiverse classes
    FitMultiverse(const FitMultiverse&&);

    std::string Checksum() const;

    std::vector<FitUniverse> fUnivs;
  };

}

// TODO split into two headers here

#include "sbnana/CAFAna/Core/ISyst.h"

namespace ana
{
  /// \brief Equivalent to FitMultiverse, but storing experiment-specific
  /// systematic type SystT
  ///
  /// Internal implementation details: we store the systematics as FitSyst* in
  /// the base class, but since we know they can only be set by the
  /// constructors, which we control, it is safe to cast them back to the
  /// derived class on-demand.
  template<class SystT> class _Multiverse: public FitMultiverse
  {
  public:
    static_assert(std::is_base_of_v<IFitSyst, SystT>);

    /// See FitMultiverse::Hypercross
    static const _Multiverse& Hypercross(const std::vector<const SystT*>& systs,
                                         int nsigma = 3)
    {
      return Cached(FitMultiverse::Hypercross(ConvertSysts(systs), nsigma));
    }

    using FitMultiverse::kTrulyRandom;

    /// See FitMultiverse::RandomGas
    static const _Multiverse& RandomGas(const std::vector<const SystT*>& systs,
                                        int Nuniv,
                                        unsigned int seed)
    {
      return Cached(FitMultiverse::RandomGas(ConvertSysts(systs), Nuniv, seed));
    }

    /// Details of a particular universe
    const _Universe<SystT>& GetUniverse(int i) const
    {
      // TODO I think this cast is safe because the layout is identical?
      return *((const _Universe<SystT>)&fUnivs[i]);
    }

  protected:
    /// Helper constructor for \ref Cached
    _Multiverse(const FitMultiverse&& m) : FitMultiverse(std::move(m)) {}

    /// Helper function for named constructors
    static std::vector<const IFitSyst*> ConvertSysts(const std::vector<const SystT*> systs)
    {
      return {systs.begin(), systs.end()};
    }

    /// \brief We need to retain the property that identically-defined
    /// multiverses are the same object (get the same pointers)
    ///
    /// If there is already a _Multiverse corresponding to \a m return that,
    /// otherwise construct one
    static const _Multiverse& Cached(const FitMultiverse& m)
    {
      // These will be deleted by the underlying FitMultiverse, so no need for
      // unique_ptr here
      static std::unordered_map<const FitMultiverse*, const _Multiverse*> cache;
      if(!cache.count(&m)){
        // We need to create a _Multiverse inheriting from the passed-in
        // FitMultiverse. Trouble is, we mustn't create a new FitMultiverse
        // with the same definition, so we have to move the input to become the
        // base of this new object. This is OK because 1. we only just created
        // this FitMultiverse in the named constructor and 2. FitMultiverse
        // defines a move constructor for us that takes care of Registry
        // registration etc.
        cache.emplace(&m, new _Multiverse(std::move(m)));
      }
      return *cache[&m];
    }
  };

  using Universe = _Universe<ISyst>;
  using Multiverse = _Multiverse<ISyst>;
}
