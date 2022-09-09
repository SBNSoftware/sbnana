#pragma once

#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"

#include "sbnana/CAFAna/Core/IRecordSource.h"
#include "sbnana/CAFAna/Core/IRecordSink.h"

#include "cafanacore/Passthrough.h"

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
  class SpectrumLoader: public SpectrumLoaderBase, public Passthrough<caf::SRSpillProxy>
  {
  public:
    SpectrumLoader(const std::string& wildcard, int max = 0);
    SpectrumLoader(const std::vector<std::string>& fnames, int max = 0);
    /// Named constructor for SAM projects
    static SpectrumLoader* FromSAMProject(const std::string& proj,
                                          int fileLimit = -1);
    virtual ~SpectrumLoader();

    virtual void Go() override;

    virtual void PrintGraph(std::ostream& os) const override;
    // Print to stdout
    virtual void PrintGraph() const;

  protected:
    SpectrumLoader();

    // Move operations
    SpectrumLoader(SpectrumLoader&&) = default;
    SpectrumLoader& operator=(SpectrumLoader&&) = default;

    // No copy operations because I don't want to deal with pointers
    SpectrumLoader(const SpectrumLoader&) = delete;
    SpectrumLoader& operator=(const SpectrumLoader&) = delete;

    virtual void HandleFile(TFile* f, Progress* prog = 0);

    /// Save accumulated exposures into the individual spectra
    virtual void StoreExposures();

    int max_entries;
  };


  class NullLoader: public SpectrumLoader
  {
  public:
    virtual void Go() override {}
  };
  static NullLoader kNullLoader;


  static NullSource<caf::SRSpillProxy> kNullSpillSource;
  static NullSource<caf::SRSliceProxy> kNullSliceSource;
}
