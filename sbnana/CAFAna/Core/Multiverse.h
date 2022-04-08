#pragma once

#include "CAFAna/Core/INamed.h"

#include "sbnana/CAFAna/Core/SystShifts.h"

class TDirectory;

namespace ana
{
  /// \brief Collection of "universes" (SystShifts)
  ///
  /// Two Multiverses are equivalent if-and-only-if they have the same address
  ///
  /// Universe zero is always the nominal universe (no systematic shifts)
  class Multiverse: public INamed
  {
  public:
    /// \brief Named constructor. Scan each parameter sequentially
    ///
    /// \param systs  The list of systematic parameters to scan
    /// \param nsigma Each parameter will be scanned from -nsigma to +nsigma
    static const Multiverse& Hypercross(const std::vector<const ISyst*>& systs,
                                        int nsigma = 3);

    static const unsigned int kTrulyRandom = 123456789;
    
    /// \brief Named constructor. Throw all parameters as gaussians
    ///
    /// \param systs The list of systematic parameters to vary
    /// \param Nuniv Number of universes to generate
    /// \param seed  Pseudo-random seed. Pass Multiverse::kTrulyRandom for a
    ///              real random seed. Beware incompatibility between universes
    ///              with different seeds.
    static const Multiverse& RandomGas(const std::vector<const ISyst*>& systs,
                                       int Nuniv,
                                       unsigned int seed);

    /// Total number of universes, including nominal at index 0
    unsigned int NUniv() const {return fUnivs.size();}

    /// Details of a particular universe
    const SystShifts& Universe(int i) const {return fUnivs[i];}

    void SaveTo(TDirectory* dir, const std::string& name) const;

    // Usually these return unique_ptr, but Multiverses are globally managed
    static const Multiverse* LoadFrom(TDirectory* dir, const std::string& name);

  protected:
    Multiverse(const std::string& shortName,
               const std::string& latexName,
               const std::vector<SystShifts>& univs);

    Multiverse(const Multiverse&) = delete;
    Multiverse(const Multiverse&&) = delete;

    std::string Checksum() const;

    // TODO disentangle SystShifts vs SystParams
    std::vector<SystShifts> fUnivs;
  };
}
