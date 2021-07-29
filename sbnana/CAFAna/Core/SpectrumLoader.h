#pragma once

#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"

#include "sbnana/CAFAna/Core/IRecordSource.h"

class TFile;

namespace ana
{
  class Progress;

  /// \brief Collaborates with \ref Spectrum and \ref OscillatableSpectrum to
  /// fill spectra from CAF files.
  ///
  /// The outside user should just pass a filename, wildcard, SAM dataset name,
  /// or query to the constructor. Then construct all the spectra you
  /// need. They will register with this loader. Finally, calling \ref Go will
  /// cause all the spectra to be filled at once. After this the loader may not
  /// be used again.
  class SpectrumLoader: public SpectrumLoaderBase, public ISpillSource, public ISliceSource
  {
  public:
    SpectrumLoader(const std::string& wildcard, int max = 0);
    SpectrumLoader(const std::vector<std::string>& fnames, int max = 0);
    /// Named constructor for SAM projects
    static SpectrumLoader* FromSAMProject(const std::string& proj,
                                          int fileLimit = -1);
    virtual ~SpectrumLoader();

    virtual void Go() override;

    using ISpillSource::GetVar;
    using ISliceSource::GetVar;
    using ISpillSource::operator[];
    using ISliceSource::operator[];

  protected:
    SpectrumLoader();

    // Move operations
    SpectrumLoader(SpectrumLoader&&) = default;
    SpectrumLoader& operator=(SpectrumLoader&&) = default;

    // No copy operations because I don't want to deal with pointers
    SpectrumLoader(const SpectrumLoader&) = delete;
    SpectrumLoader& operator=(const SpectrumLoader&) = delete;

    virtual void HandleFile(TFile* f, Progress* prog = 0);

    virtual void HandleRecord(caf::SRSpillProxy* sr);

    /// Save accumulated exposures into the individual spectra
    virtual void StoreExposures();

    /// All unique cuts contained in fHistDefs
    //    std::vector<Cut> fAllCuts;
    //    std::vector<double> fLivetimeByCut; ///< Indexing matches fAllCuts
    //    std::vector<double> fPOTByCut;      ///< Indexing matches fAllCuts
    int max_entries;
  };

  // TODO does this actually work right?
  /// \brief Dummy loader that doesn't load any files
  ///
  /// Useful when a loader is required for a component you want to ignore
  class NullLoader: public SpectrumLoader
  {
  public:
    NullLoader();
    virtual ~NullLoader();
    virtual void Go() override;
  };
  /// \brief Dummy loader that doesn't load any files
  ///
  /// Useful when a loader is required for a component you want to ignore
  static NullLoader kNullLoader;
}
