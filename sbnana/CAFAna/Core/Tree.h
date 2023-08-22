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
    virtual void SaveTo( TDirectory* dir ) const;
  protected:
    std::map< std::string, std::vector<double>> fBranchEntries;
    std::string fTreeName;
    std::vector<std::string> fOrderedBranchNames;
    long long fNEntries;
    double fPOT;
    double fLivetime;
    bool fSaveRunSubEvt;
    bool fSaveSliceNum;
  };

  // Similar to Tree but for event weights e.g. to make splines...
  class WeightsTree
  {
  public:
    WeightsTree( const std::string name, const std::vector<std::string>& labels,
                 const unsigned int nSigma, const bool saveRunSubEvt, const bool saveSliceNum, const unsigned int nWeightsExpected );
    /// Function to update protected members (the branches). DO NOT USE outside of the filling.
    void UpdateEntries ( const std::map<std::string, std::vector<double>> valsMap, const std::map<std::string, std::vector<double>> weightMap );
    /// Function to update protected members (the exposures). DO NOT USE outside of the filling.
    void UpdateExposure ( const double pot, const double livetime );
    // Utilities
    double POT() const {return fPOT;} // as in Spectrum
    double Livetime() const {return fLivetime;} // as in Spectrum
    long long Entries() const {return fNEntries;}
    int NSigma() const {return fNSigma;} // return as an int, then -NSigma to NSigma means something
    bool SaveRunSubEvent() const {return fSaveRunSubEvt;}
    bool SaveSliceNum() const {return fSaveSliceNum;}
    void OverridePOT(double newpot) {fPOT = newpot;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    void OverrideLivetime(double newlive) {fLivetime = newlive;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    virtual void SaveTo( TDirectory* dir ) const = 0;
  protected:
    std::map< std::string, std::vector<double>> fBranchEntries;
    std::map< std::string, std::vector<std::vector<double>>> fBranchWeightEntries;
    std::string fTreeName;
    std::vector<std::string> fOrderedBranchNames;
    std::vector<std::string> fOrderedBranchWeightNames;
    long long fNEntries;
    double fPOT;
    double fLivetime;
    bool fSaveRunSubEvt;
    bool fSaveSliceNum;
    unsigned int fNSigma;
    unsigned int fNWeightsExpected;
  };

  class NSigmasTree : public WeightsTree
  {
  public:
    /// constructor with a vector of \ref ISyst
    NSigmasTree( const std::string name, const std::vector<std::string>& labels,
                 SpectrumLoaderBase& loader,
                 const std::vector<const ISyst*>& systsToStore, const SpillCut& spillcut,
                 const Cut& cut, const SystShifts& shift = kNoShift, const unsigned int nSigma = 3, const bool saveRunSubEvt = false, const bool saveSliceNum = false );
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
                    const std::vector<std::vector<Var>>& univsKnobs, const unsigned int nUniverses,
                    const SpillCut& spillcut,
                    const Cut& cut, const SystShifts& shift = kNoShift, const bool saveRunSubEvt = false, const bool saveSliceNum = false );
    void SaveTo( TDirectory* dir ) const override;
    unsigned int NUniverses() const {return fNSigma;}
  };

}
