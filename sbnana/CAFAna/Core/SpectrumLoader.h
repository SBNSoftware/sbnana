#pragma once

#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"

#include "sbnana/CAFAna/Core/IRecordSource.h"
#include "sbnana/CAFAna/Core/IRecordSink.h"

class TFile;

namespace ana
{
  class Progress;

  // TODO drop the "SBN" and define all this in RecordSource?
  /// SBNSpillSource is special since it also knows how to loop over slices
  class SBNSpillSource: public ISpillSink, public ISpillSource, public ISliceSource
  {
  public:

    using ISliceSource::GetVar;
    using ISliceSource::GetVars;
    beta::IValueSource& operator[](const Var& var){return GetVar(var);}

    using ISliceSource::GetCut;
    ISliceSource& operator[](const Cut& cut){return ISliceSource::GetCut(cut);}

    SBNSpillSource& GetCut(const SpillCut& cut);

    SBNSpillSource& operator[](const SpillCut& cut){return GetCut(cut);}

    using ISpillSource::GetVar;
    using ISpillSource::GetVars;
    beta::IValueSource& operator[](const SpillVar& var){return GetVar(var);}

    ISliceSource& Ensemble(const std::vector<Weight>& weis, int multiverseId)
    {
      return ISliceSource::Ensemble(weis, multiverseId);
    }

    ISpillSource& Ensemble(const std::vector<SpillWeight>& weis, int multiverseId)
    {
      return ISpillSource::Ensemble(weis, multiverseId);
    }

    virtual void HandleRecord(const caf::SRSpillProxy* spill, double weight, int universeId) override;
    virtual void HandleEnsemble(const caf::SRSpillProxy* spill, const std::vector<double>& weights, int multiverseId) override;

    virtual void HandlePOT(double pot) override;

    virtual void HandleLivetime(double livetime) override;

    virtual unsigned int NSinks() const override;

  protected:
    SBNSpillSource(){}
  };

  /// \brief Collaborates with \ref Spectrum and \ref OscillatableSpectrum to
  /// fill spectra from CAF files.
  ///
  /// The outside user should just pass a filename, wildcard, SAM dataset name,
  /// or query to the constructor. Then construct all the spectra you
  /// need. They will register with this loader. Finally, calling \ref Go will
  /// cause all the spectra to be filled at once. After this the loader may not
  /// be used again.
  class SpectrumLoader: public SpectrumLoaderBase, public SBNSpillSource
  {
  public:
    SpectrumLoader(const std::string& wildcard, int max = 0);
    SpectrumLoader(const std::vector<std::string>& fnames, int max = 0);
    /// Named constructor for SAM projects
    static SpectrumLoader* FromSAMProject(const std::string& proj,
                                          int fileLimit = -1);
    virtual ~SpectrumLoader();

    virtual void Go() override;

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
