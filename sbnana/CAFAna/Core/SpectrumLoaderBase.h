#pragma once

#include <cassert>
#include <functional>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "cafanacore/IFileSource.h"
#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Weight.h"

namespace caf{class StandardRecord;}

class TFile;

namespace ana
{
  class Spectrum;
  class ReweightableSpectrum;

  /// Base class for the various types of spectrum loader
  class SpectrumLoaderBase
  {
  public:
    virtual ~SpectrumLoaderBase();

    /// Load all the registered spectra
    virtual void Go() = 0;

    /// Indicate whether or not \ref Go has been called
    virtual bool Gone() const {return fGone;}

  protected:
    /// Component of other constructors
    SpectrumLoaderBase();
    /// Construct from a filename, wildcard, SAM definition, or SAM query
    SpectrumLoaderBase(const std::string& wildcard);
    /// Construct from an explicit list of files
    SpectrumLoaderBase(const std::vector<std::string>& fnames);

    // Move operations
    SpectrumLoaderBase(SpectrumLoaderBase&&) = default;
    SpectrumLoaderBase& operator=(SpectrumLoaderBase&&) = default;

    // No copy operations because I don't want to deal with pointers
    SpectrumLoaderBase(const SpectrumLoaderBase&) = delete;
    SpectrumLoaderBase& operator=(const SpectrumLoaderBase&) = delete;

    /// Figure out if \a str is a wildcard or SAM query and return a source
    IFileSource* WildcardOrSAMQuery(const std::string& str) const;

    /// Forwards to \ref fFileSource
    int NFiles() const;

    TFile* GetNextFile();

    std::string fWildcard;
    std::unique_ptr<IFileSource> fFileSource;

    bool fGone; ///< Has Go() been called? Can't add more histograms after that

    double fPOT;
    double fPOTFromHist; ///< Accumulated by calls to \ref GetNextFile
    int fNReadouts;
  };
}
