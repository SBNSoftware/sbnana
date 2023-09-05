#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Cut.h"

#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include <memory>
#include <string>
#include <vector>

class TDirectory;
class TTree;
class TBranch;

namespace ana
{

  class Tree
  {
  public:
    /// constructor with a vector of \ref Var
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<Var>& vars, const SpillCut& spillcut,
          const Cut& cut, const SystShifts& shift = kNoShift, const bool saveRunSubEvt = false, const bool saveSliceNum = false );
    /// constructor with a vector of \ref MultiVar
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<MultiVar>& vars, const SpillCut& spillcut,
          const Cut& cut, const SystShifts& shift = kNoShift, const bool saveRunSubEvt = false, const bool saveSliceNum = false );
    /// constructor with a vector of \ref SpillVar
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<SpillVar>& vars, const SpillCut& spillcut, const bool saveRunSubEvt = false );
    /// constructor with a vector of \ref SpillMultiVar
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<SpillMultiVar>& vars, const SpillCut& spillcut, const bool saveRunSubEvt = false );
    /// constructor with a vector of \ref TruthVar
    /// Dedicated for signal efficiency, thus has slightly different behavior
    /// - spillcut: It does not cut event by spillcut, but creates a branch "SpillCutType"
    /// - truthcut: Signal definition is defined here. Only the TrueInteraction (=nu) that satisfies truthcut will be filled
    /// - SignalSelection: Reco cut that defined your signal selection. 
    ///   For each nu, it loops over reco slices, and find if any matched slice pass this cut, and save this to "CutType" branch
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<TruthVar>& vars, const SpillCut& spillcut,
          const TruthCut& truthcut,
          const Cut& SignalSelection,
          const SystShifts& shift = kNoShift,
          const bool saveRunSubEvt = false );

    // Add functionality to update the protected stuff from elsewhere
    /// Function to update protected members (the branches). DO NOT USE outside of the filling.
    virtual void UpdateEntries ( const std::map<std::string, std::vector<double>> valsMap );
    /// Function to update protected members (the exposures). DO NOT USE outside of the filling.
    void UpdateExposure ( const double pot, const double livetime );
    // Utilities
    double POT() const {return fPOT;} // as in Spectrum
    double Livetime() const {return fLivetime;} // as in Spectrum
    long long Entries() const {return fNEntries;}
    bool SaveRunSubEvent() const {return fSaveRunSubEvt;}
    bool SaveSliceNum() const {return fSaveSliceNum;}
    void OverridePOT(double newpot) {fPOT = newpot;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    void OverrideLivetime(double newlive) {fLivetime = newlive;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    bool SaveTruthCutType() const {return fSaveTruthCutType;}
    Cut GetSignalSelectionCut() const {return SignalSelection;}
    virtual void SaveTo( TDirectory* dir ) const;
  protected:
    friend class WeightsTree;
    std::map< std::string, std::vector<double>> fBranchEntries;
    std::string fTreeName;
    std::vector<std::string> fOrderedBranchNames;
    long long fNEntries;
    double fPOT;
    double fLivetime;
    bool fSaveRunSubEvt;
    bool fSaveSliceNum;
    bool fSaveTruthCutType;
    const Cut SignalSelection;
  };

  /// Allows to make covariance matrices using \ref Tree objects cleverly
  class CovarianceMatrixTree : public Tree
  {
  public:
    /// constructor with a vector \ref SpillMultiVar corresponding to the items to be used in binning decision and the weights
    CovarianceMatrixTree( const std::string name, const std::vector<std::string>& labels, const std::string weightLabel,
                          const std::vector< std::vector<std::pair<std::string,std::pair<double,double>>> >& binning,
                          SpectrumLoaderBase& loader, const std::vector<SpillMultiVar>& vars, const unsigned int nUniverses,
                          const SpillCut& spillcut );
    /// Function to update protected members (the branches). DO NOT USE outside of the filling.
    void UpdateEntries ( const std::map<std::string, std::vector<double>> valsMap ) override;
    void SaveTo( TDirectory* dir ) const override;
    void PrintBinning() const;
  private:
    std::map< unsigned int, std::vector< std::pair<std::string,std::pair<double,double>> >> fBinning;
    std::vector<std::vector<double>> fWeights;
    std::string fWeightLabel;
    unsigned int fNUniverses;
  };

  // Similar to Tree but for event weights e.g. to make splines...
  class WeightsTree
  {
  public:
    WeightsTree( const std::string name, const std::vector<std::string>& labels,
                 const std::vector<unsigned int>& nWeights, const bool saveRunSubEvt, const bool saveSliceNum );
    WeightsTree( const std::string name, const std::vector<std::string>& labels,
                 const std::vector<std::pair<int,int>>& nSigma, const bool saveRunSubEvt, const bool saveSliceNum );
    /// Function to merge the WeightsTree with a Tree, assuming both have fSaveRunSubEvt and fSaveSliceNum set to true
    void MergeTree( const Tree& inTree );
    /// Function to update protected members (the branches). DO NOT USE outside of the filling.
    void UpdateEntries ( const std::map<std::string, std::vector<double>> valsMap, const std::map<std::string, std::vector<double>> weightMap );
    /// Function to update protected members (the exposures). DO NOT USE outside of the filling.
    void UpdateExposure ( const double pot, const double livetime );
    // Utilities
    double POT() const {return fPOT;} // as in Spectrum
    double Livetime() const {return fLivetime;} // as in Spectrum
    long long Entries() const {return fNEntries;}
    int NSigmaLo(const std::string& systToNSigma) const {return fNSigmasLo.at(systToNSigma);}
    int NSigmaHi(const std::string& systToNSigma) const {return fNSigmasHi.at(systToNSigma);}
    unsigned int NWeights(const std::string& systToNWts) const {return fNWeightsExpected.at(systToNWts);}
    bool SaveRunSubEvent() const {return fSaveRunSubEvt;}
    bool SaveSliceNum() const {return fSaveSliceNum;}
    void OverridePOT(double newpot) {fPOT = newpot;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    void OverrideLivetime(double newlive) {fLivetime = newlive;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    virtual void SaveTo( TDirectory* dir ) const = 0;
  protected:
    std::map< std::string, std::vector<double>> fBranchEntries;
    std::map< std::string, std::vector<std::vector<double>>> fBranchWeightEntries;
    std::map< std::string, unsigned int> fNSigmasLo;
    std::map< std::string, unsigned int> fNSigmasHi;
    std::map< std::string, unsigned int> fNWeightsExpected;
    //std::map< std::string, std::vector<double>> fNWeightsThrows; // stores the expected ordering of NSigmas per knob or universe throws per knob.
    std::string fTreeName;
    std::vector<std::string> fOrderedBranchNames;
    std::vector<std::string> fOrderedBranchWeightNames;
    long long fNEntries;
    double fPOT;
    double fLivetime;
    bool fSaveRunSubEvt;
    bool fSaveSliceNum;
  };

  class NSigmasTree : public WeightsTree
  {
  public:
    /// constructor with a vector of \ref ISyst
    NSigmasTree( const std::string name, const std::vector<std::string>& labels,
                 SpectrumLoaderBase& loader,
                 const std::vector<const ISyst*>& systsToStore, const std::vector<std::pair<int,int>>& nSigma,
                 const SpillCut& spillcut,
                 const Cut& cut, const SystShifts& shift = kNoShift, const bool saveRunSubEvt = false, const bool saveSliceNum = false );
    /// constructor with a vector of \ref ISyst, but for TrueTree
    NSigmasTree( const std::string name, const std::vector<std::string>& labels,
                 SpectrumLoaderBase& loader,
                 const std::vector<const ISyst*>& systsToStore, const std::vector<std::pair<int,int>>& nSigma,
                 const TruthCut& truthcut, const SystShifts& shift = kNoShift, const bool saveRunSubEvt = false);
    void SaveTo( TDirectory* dir ) const override;
    void SaveToSplines( TDirectory* dir ) const;
    void SaveToGraphs( TDirectory* dir ) const;
  };

  class NUniversesTree : public WeightsTree
  {
  public:
    /// constructor with a vector of vectors of \ref Var corresponding to a number of universes for which we want to extract weights
    NUniversesTree( const std::string name, const std::vector<std::string>& labels,
                    SpectrumLoaderBase& loader,
                    const std::vector<std::vector<Var>>& univsKnobs, const std::vector<unsigned int>& nUniverses,
                    const SpillCut& spillcut,
                    const Cut& cut, const SystShifts& shift = kNoShift, const bool saveRunSubEvt = false, const bool saveSliceNum = false );
    void SaveTo( TDirectory* dir ) const override;
  };

}
