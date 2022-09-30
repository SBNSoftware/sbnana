#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"

#include <set>

namespace ana
{
  class Progress;

  // Example reduction steps
  void ClearTrueParticles(caf::StandardRecord* sr);

  /// \brief Create smaller CAFs
  ///
  /// This class produces new CAFs, removing entries that fail a cut. It also
  /// allows the event record to be edited in custom ways.
  class FileReducer: protected SpectrumLoaderBase
  {
  public:
    FileReducer(const std::string& wildcard,
		const std::string& outfile);
    FileReducer(const std::vector<std::string>& fnames,
                const std::string& outfile);
    virtual ~FileReducer();

    /// Only copy records to the output file if they pass this cut
    void AddSpillCut(const SpillCut& cut);
    void AddSliceCut(const SliceCut& cut);

    /// \brief If called, only events whose run/subrun/event occur in \a fname
    /// will be retained.
    void SetEventList(const std::string& fname);

    typedef void (ReductionFunc)(caf::StandardRecord*);
    //    typedef void (ReductionFuncWithProxy)(caf::StandardRecord*,
    //                                          const caf::Proxy<caf::StandardRecord>*);

    /// Run the specified reduction function over each event
    void AddReductionStep(const std::function<ReductionFunc> &f) {fReductionFuncs.push_back(f);}
    //    void AddReductionStep(const std::function<ReductionFuncWithProxy> &f) {fReductionFuncsWithProxy.push_back(f);}

    /// Override any metadata key in the output file
    void SetMetadata(const std::string& key, const std::string& val)
    {
      // Let's assume the user is just setting string keys
      fMetaMap[key] = "\""+val+"\"";
    }

    virtual void Go() override;

  protected:
    void HandleFile(TFile* fin, TFile* fout, TTree*& trOut, Progress* prog,
                    long& nRecSeen, long& nRecPassed);

    void HandleNestedTree(TFile* fout, TTree* recTree, TTree*& trOut,
                          Progress* prog,
                          long& nRecSeen, long& nRecPassed);

    void HandleFlatTree(TFile* fout, TTree* recTree, TTree*& trOut,
                        Progress* prog,
                        long& nRecSeen, long& nRecPassed);

    void CopyGlobalTree(TFile* fin, TFile* fout);

    void UpdateMetadata(std::map<std::string, std::string>& meta,
                        const std::set<std::string>& mask,
                        const std::vector<std::string>& fnames) const;

    /// Strip all information out of this record and tag it as a husk
    void Huskify(caf::StandardRecord* sr) const;

    std::string fOutfile;
    SpillCut*   fSpillCut;
    SliceCut*   fSliceCut;

    std::set<std::tuple<int, int, int>> fEventList;

    std::vector<std::function<ReductionFunc>> fReductionFuncs;
    //    std::vector<std::function<ReductionFuncWithProxy>> fReductionFuncsWithProxy;

    std::map<std::string, std::string> fMetaMap;
  };
}
