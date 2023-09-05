#pragma once

#include <functional>
#include <list>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/IFileSource.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Core/Var.h"

namespace caf{class StandardRecord;}

class TFile;
class TH1;
class TTree;

namespace ana
{
  class Spectrum;
  class ReweightableSpectrum;
  class Tree;
  class WeightsTree;
  class NSigmasTree;
  class NUniversesTree;

  /// Is this data-file representing beam spills or cosmic spills?
  enum DataSource{
    // TODO there are no longer any actions taken as a result of this
    // distinction. Remove after SA.
    kBeam,
    kCosmic
  };

  /// Base class for the various types of spectrum loader
  class SpectrumLoaderBase
  {
  public:

    friend class ReweightableSpectrum;
    friend class NDOscillatableSpectrum;
    friend class OscillatableSpectrum;
    friend class Spectrum;

    virtual ~SpectrumLoaderBase();

    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const Var& var,
                             const SpillCut& spillcut,
                             const Cut& cut,
                             const SystShifts& shift,
                             const Var& wei = kUnweighted);

    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const MultiVar& var,
                             const SpillCut& spillcut,
                             const Cut& cut,
                             const SystShifts& shift,
                             const Var& wei = kUnweighted);

    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const SpillVar& var,
                             const SpillCut& cut,
                             const SpillVar& wei = kSpillUnweighted);

    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const SpillMultiVar& var,
                             const SpillCut& cut,
                             const SpillVar& wei = kSpillUnweighted);

    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const TruthVar& var,
                             const TruthCut truthcut,
                             const SpillCut& spillcut,
                             const SystShifts& shift,
                             const TruthVar& wei = kTruthUnweighted);
    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const TruthMultiVar& var,
                             const TruthCut truthcut,
                             const SpillCut& spillcut,
                             const SystShifts& shift,
                             const TruthVar& wei = kTruthUnweighted);
    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const TruthVar& var,
                             const TruthCut truthcut,
                             const SpillCut& spillcut,
                             const Cut& cut, // loop over reco slices and see if any matched to this truth and pass "cut"
                             const SystShifts& shift,
                             const TruthVar& wei = kTruthUnweighted);
    /// For use by the \ref Spectrum constructor
    virtual void AddSpectrum(Spectrum& spect,
                             const TruthMultiVar& var,
                             const TruthCut truthcut,
                             const SpillCut& spillcut,
                             const Cut& cut, // loop over reco slices and see if any matched to this truth and pass "cut"
                             const SystShifts& shift,
                             const TruthVar& wei = kTruthUnweighted);

    /// For use by the constructors of \ref ReweightableSpectrum subclasses
    virtual void AddReweightableSpectrum(ReweightableSpectrum& spect,
                                         const Var& var,
                                         const Cut& cut,
                                         const SystShifts& shift,
                                         const Var& wei);

    /// For use by the constructors of \ref ReweightableSpectrum subclasses
    virtual void AddReweightableSpectrum(ReweightableSpectrum& spect,
                                         const Var& var,
                                         const SpillCut& spillcut,
                                         const SliceCut& slicecut,
                                         const SystShifts& shift,
                                         const Var& wei);

    /// For use by the constructors of \ref Tree class
    virtual void AddTree(Tree& tree,
                         const std::vector<std::string>& labels,
                         const std::vector<Var>& vars,
                         const SpillCut& spillcut,
                         const Cut& cut,
                         const SystShifts& shift);

    /// For use by the constructors of \ref Tree class
    virtual void AddTree(Tree& tree,
                         const std::vector<std::string>& labels,
                         const std::vector<MultiVar>& vars,
                         const SpillCut& spillcut,
                         const Cut& cut,
                         const SystShifts& shift);

    /// For use by the constructors of \ref Tree class
    virtual void AddTree(Tree& tree,
                         const std::vector<std::string>& labels,
                         const std::vector<SpillVar>& vars,
                         const SpillCut& spillcut);

    /// For use by the constructors of \ref Tree class
    virtual void AddTree(Tree& tree,
                         const std::vector<std::string>& labels,
                         const std::vector<SpillMultiVar>& vars,
                         const SpillCut& spillcut);

    /// For use by the constructors of \ref Tree class
    virtual void AddTree(Tree& tree,
                         const std::vector<std::string>& labels,
                         const std::vector<TruthVar>& vars,
                         const SpillCut& spillcut,
                         const TruthCut& truthcut,
                         const SystShifts& shift);

    /// For use by the constructors of \ref NSigmasTree class
    virtual void AddNSigmasTree(NSigmasTree& tree,
                                const std::vector<std::string>& labels,
                                const std::vector<const ISyst*>& systsToStore,
                                const SpillCut& spillcut,
                                const Cut& cut,
                                const SystShifts& shift);

    /// For use by the constructors of \ref NUniversesTree class
    virtual void AddNUniversesTree(NUniversesTree& tree,
                                   const std::vector<std::string>& labels,
                                   const std::vector<std::vector<Var>>& univKnobs,
                                   const SpillCut& spillcut,
                                   const Cut& cut,
                                   const SystShifts& shift);

    /// Load all the registered spectra
    virtual void Go() = 0;

    /// Indicate whether or not \ref Go has been called
    virtual bool Gone() const {return fGone;}

  protected:
    /// Component of other constructors
    SpectrumLoaderBase(DataSource src = kBeam);
    /// Construct from a filename, wildcard, SAM definition, or SAM query
    SpectrumLoaderBase(const std::string& wildcard, DataSource src = kBeam);
    /// Construct from an explicit list of files
    SpectrumLoaderBase(const std::vector<std::string>& fnames, DataSource src = kBeam);

    // Move operations
    SpectrumLoaderBase(SpectrumLoaderBase&&) = default;
    SpectrumLoaderBase& operator=(SpectrumLoaderBase&&) = default;

    // No copy operations because I don't want to deal with pointers
    SpectrumLoaderBase(const SpectrumLoaderBase&) = delete;
    SpectrumLoaderBase& operator=(const SpectrumLoaderBase&) = delete;

    /// Figure out if \a str is a wildcard or SAM query and return a source
    IFileSource* WildcardOrSAMQuery(const std::string& str) const;

    friend class SpectrumLoaderMockData;
    virtual void RemoveSpectrum(Spectrum*);
    virtual void RemoveReweightableSpectrum(ReweightableSpectrum*);

    /// Forwards to \ref fFileSource
    int NFiles() const;

    TFile* GetNextFile();

    std::string fWildcard;
    std::unique_ptr<IFileSource> fFileSource;

    DataSource fSource;

    bool fGone; ///< Has Go() been called? Can't add more histograms after that

    double fPOT;
    double fPOTFromHist; ///< Accumulated by calls to \ref GetNextFile
    int fNReadouts;

    /// \brief Helper class for \ref SpectrumLoaderBase
    ///
    /// List of Spectrum and OscillatableSpectrum, some utility functions
    struct SpectList
    {
      void Erase(Spectrum* s);
      void Erase(ReweightableSpectrum* os);
      void RemoveLoader(SpectrumLoaderBase* l);
      size_t TotalSize() const;
      void GetSpectra(std::vector<Spectrum*>& ss);
      void GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss);

      std::vector<Spectrum*> spects;
      std::vector<ReweightableSpectrum*> rwSpects;
    };

    /// \brief Helper class for \ref SpectrumLoaderBase
    ///
    /// Functions like std::map<T, U> except it should be faster to iterate
    /// through the elements (while slower to fill) and it knows to compare Ts
    /// via their ID() function. Various methods that forward through to the
    /// \ref SpectList at the end of the chain.
    template<class T, class U> struct IDMap
    {
      U& operator[](const T& key);

      // Make class iterable. Keep inline for speed
      typedef typename std::vector<std::pair<T, U>>::iterator it_t;
      inline it_t begin(){return fElems.begin();}
      inline it_t end(){return fElems.end();}

      template<class V> void Erase(const V& v);
      void RemoveLoader(SpectrumLoaderBase* l);
      void Clear();
      size_t TotalSize();
      void GetSpectra(std::vector<Spectrum*>& ss);
      void GetReweightableSpectra(std::vector<ReweightableSpectrum*>& ss);
    protected:
      std::vector<std::pair<T, U>> fElems;
    };

    template<class T> class _VarOrMultiVar
    {
    public:
      // v could easily be a temporary, have to make a copy
      _VarOrMultiVar(const _Var<T>& v) : fVar(new _Var<T>(v)), fMultiVar(0) {}
      _VarOrMultiVar(const _MultiVar<T>& v) : fVar(0), fMultiVar(new _MultiVar<T>(v)) {}
      ~_VarOrMultiVar() {delete fVar; delete fMultiVar;}

      _VarOrMultiVar(const _VarOrMultiVar& v)
        : fVar(v.fVar ? new _Var<T>(*v.fVar) : 0),
          fMultiVar(v.fMultiVar ? new _MultiVar<T>(*v.fMultiVar) : 0)
      {
      }

      _VarOrMultiVar(_VarOrMultiVar&& v)
      {
        fVar = v.fVar;
        fMultiVar = v.fMultiVar;
        v.fVar = 0;
        v.fMultiVar = 0;
      }

      bool IsMulti() const {return fMultiVar;}
      const _Var<T>& GetVar() const {assert(fVar); return *fVar;}
      const _MultiVar<T>& GetMultiVar() const {assert(fMultiVar); return *fMultiVar;}

      int ID() const {return fVar ? fVar->ID() : fMultiVar->ID();}

    protected:
      const _Var<T>* fVar;
      const _MultiVar<T>* fMultiVar;
    };

    typedef _VarOrMultiVar<caf::SRSliceProxy> VarOrMultiVar;
    typedef _VarOrMultiVar<caf::SRSpillProxy> SpillVarOrMultiVar;
    typedef _VarOrMultiVar<caf::SRTrueInteractionProxy> TruthVarOrMultiVar;

    template<class T>
    friend bool operator<(const _VarOrMultiVar<T>& a, const _VarOrMultiVar<T>& b) noexcept;

    /// \brief All the spectra that need to be filled
    ///
    /// [spillcut][shift][cut][wei][var]
    IDMap<SpillCut, IDMap<SystShifts, IDMap<Cut, IDMap<Var, IDMap<VarOrMultiVar, SpectList>>>>> fHistDefs;
    /// [spillcut][spillwei][spillvar]
    IDMap<SpillCut, IDMap<SpillVar, IDMap<SpillVarOrMultiVar, SpectList>>> fSpillHistDefs;

    /// [spillcut][shift][truthcut][truthweight][truthvar]
    IDMap<SpillCut, IDMap<SystShifts, IDMap<TruthCut, IDMap<TruthVar, IDMap<TruthVarOrMultiVar, SpectList>>>>> fTruthHistDefs;
    /// [spillcut][cut][shift][truthcut][truthweight][truthvar]
    IDMap<SpillCut, IDMap<Cut, IDMap<SystShifts, IDMap<TruthCut, IDMap<TruthVar, IDMap<TruthVarOrMultiVar, SpectList>>>>>> fTruthHistWithCutDefs;

    // TODO: Probably someone can make a more efficient version of SpectList
    //       that works with Tree objects... In the meantime, let's use a standard
    //       map. But otherwise, let's keep it the same way...
    std::map<SpillCut, std::map<SystShifts, std::map<Cut, std::map<Tree*, std::map<VarOrMultiVar, std::string>>>>> fTreeDefs;
    std::map<SpillCut, std::map<Tree*, std::map<SpillVarOrMultiVar, std::string>>> fSpillTreeDefs;
    std::map<SpillCut, std::map<SystShifts, std::map<TruthCut, std::map<Tree*, std::map<TruthVarOrMultiVar, std::string>>>>> fTruthTreeDefs;
    // And a version that saves up the syst weights used to make event-by-event splines
    std::map<SpillCut, std::map<SystShifts, std::map<Cut, std::map<NSigmasTree*, std::map<const ISyst*, std::string>>>>> fNSigmasTreeDefs;
    // And a version that saves up universe-based systematic weights to make event-by-event weight lists
    std::map<SpillCut, std::map<SystShifts, std::map<Cut, std::map<NUniversesTree*, std::map<std::vector<VarOrMultiVar>, std::string>>>>> fNUniversesTreeDefs;

  };

  // For map-making
  template<class T> bool operator<(const SpectrumLoaderBase::_VarOrMultiVar<T>& a, const SpectrumLoaderBase::_VarOrMultiVar<T>& b) noexcept {return a.ID() < b.ID();}

  /// \brief Dummy loader that doesn't load any files
  ///
  /// Useful when a loader is required for a component you want to ignore
  class NullLoader: public SpectrumLoaderBase
  {
  public:
    NullLoader() {}
    ~NullLoader();

    virtual void Go() override;
    void AddSpectrum(Spectrum& spect,
                     const Var& var,
                     const SpillCut& spillcut,
                     const Cut& cut,
                     const SystShifts& shift,
                     const Var& wei = kUnweighted) override {}
    void AddSpectrum(Spectrum& spect,
                     const MultiVar& var,
                     const SpillCut& spillcut,
                     const Cut& cut,
                     const SystShifts& shift,
                     const Var& wei = kUnweighted) override {}

    void AddSpectrum(Spectrum& spect,
                     const SpillVar& var,
                     const SpillCut& cut,
                     const SpillVar& wei = kSpillUnweighted) override {}

    void AddSpectrum(Spectrum& spect,
                     const SpillMultiVar& var,
                     const SpillCut& cut,
                     const SpillVar& wei = kSpillUnweighted) override {}

    void AddSpectrum(Spectrum& spect,
                     const TruthVar& var,
                     const TruthCut truthcut,
                     const SpillCut& spillcut,
                     const SystShifts& shift,
                     const TruthVar& wei = kTruthUnweighted) override {}

    void AddSpectrum(Spectrum& spect,
                     const TruthMultiVar& var,
                     const TruthCut truthcut,
                     const SpillCut& spillcut,
                     const SystShifts& shift,
                     const TruthVar& wei = kTruthUnweighted) override {}

    void AddSpectrum(Spectrum& spect,
                     const TruthVar& var,
                     const TruthCut truthcut,
                     const SpillCut& spillcut,
                     const Cut& cut, // loop over reco slices and see if any matched to this truth and pass "cut"
                     const SystShifts& shift,
                     const TruthVar& wei = kTruthUnweighted) override {}

    void AddSpectrum(Spectrum& spect,
                     const TruthMultiVar& var,
                     const TruthCut truthcut,
                     const SpillCut& spillcut,
                     const Cut& cut, // loop over reco slices and see if any matched to this truth and pass "cut"
                     const SystShifts& shift,
                     const TruthVar& wei = kTruthUnweighted) override {}



    void AddReweightableSpectrum(ReweightableSpectrum& spect,
                                 const Var& var,
                                 const Cut& cut,
                                 const SystShifts& shift,
                                 const Var& wei) override {}

    void AddReweightableSpectrum(ReweightableSpectrum& spect,
                                 const Var& var,
                                 const SpillCut& spillcut,
                                 const Cut& cut,
                                 const SystShifts& shift,
                                 const Var& wei) override {}

    void AddTree(Tree& tree,
                 const std::vector<std::string>& labels,
                 const std::vector<Var>& vars,
                 const SpillCut& spillcut,
                 const Cut& cut,
                 const SystShifts& shift) override {}

    void AddTree(Tree& tree,
                 const std::vector<std::string>& labels,
                 const std::vector<MultiVar>& vars,
                 const SpillCut& spillcut,
                 const Cut& cut,
                 const SystShifts& shift) override {}

    void AddTree(Tree& tree,
                 const std::vector<std::string>& labels,
                 const std::vector<SpillVar>& vars,
                 const SpillCut& spillcut) override {}

    void AddTree(Tree& tree,
                 const std::vector<std::string>& labels,
                 const std::vector<SpillMultiVar>& vars,
                 const SpillCut& spillcut) override {}

    void AddTree(Tree& tree,
                 const std::vector<std::string>& labels,
                 const std::vector<TruthVar>& vars,
                 const SpillCut& spillcut,
                 const TruthCut& truthcut,
                 const SystShifts& shift) override {}

    void AddNSigmasTree(NSigmasTree& tree,
                        const std::vector<std::string>& labels,
                        const std::vector<const ISyst*>& systsToStore,
                        const SpillCut& spillcut,
                        const Cut& cut,
                        const SystShifts& shift) override {}

    void AddNUniversesTree(NUniversesTree& tree,
                           const std::vector<std::string>& labels,
                           const std::vector<std::vector<Var>>& univKnobs,
                           const SpillCut& spillcut,
                           const Cut& cut,
                           const SystShifts& shift) override {}
  };
  /// \brief Dummy loader that doesn't load any files
  ///
  /// Useful when a loader is required for a component you want to ignore
  static NullLoader kNullLoader;
}
