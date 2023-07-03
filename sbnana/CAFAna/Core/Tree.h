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
          const Cut& cut, const SystShifts& shift = kNoShift );
    /// constructor with a vector of \ref MultiVar
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<MultiVar>& vars, const SpillCut& spillcut,
          const Cut& cut, const SystShifts& shift = kNoShift );
    /// constructor with a vector of \ref SpillVar
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<SpillVar>& vars, const SpillCut& spillcut );
    /// constructor with a vector of \ref SpillMultiVar
    Tree( const std::string name, const std::vector<std::string>& labels,
          SpectrumLoaderBase& loader,
          const std::vector<SpillMultiVar>& vars, const SpillCut& spillcut );
    // Add functionality to update the protected stuff from elsewhere
    /// Function to update protected members (the branches). DO NOT USE outside of the filling.
    void UpdateEntries ( const std::map<std::string, std::vector<double>> valsMap );
    /// Function to update protected members (the exposures). DO NOT USE outside of the filling.
    void UpdateExposure ( const double pot, const double livetime );
    // Utilities
    double POT() const {return fPOT;} // as in Spectrum
    double Livetime() const {return fLivetime;} // as in Spectrum
    long long Entries() const {return fNEntries;}
    void OverridePOT(double newpot) {fPOT = newpot;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    void OverrideLivetime(double newlive) {fLivetime = newlive;} // as in Spectrum: DO NOT USE UNLESS CERTAIN THERE ISN'T A BETTER WAY!
    void SaveTo( TDirectory* dir ) const;
  protected:
    std::map< std::string, std::vector<double>> fBranchEntries;
    std::string fTreeName;
    std::vector<std::string> fOrderedBranchNames;
    long long fNEntries;
    double fPOT;
    double fLivetime;
  };

}