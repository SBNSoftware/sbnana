#include "sbnana/CAFAna/Core/Tree.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1D.h"

#include <cassert>
#include <cmath>

namespace ana
{

  //----------------------------------------------------------------------
  // Constructor for a set of Vars, typical usage for a Selected Tree
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<Var>& vars, const SpillCut& spillcut,
              const Cut& cut, const SystShifts& shift, const bool saveRunSubEvt )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt)
  {
    assert( labels.size() == vars.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      fOrderedBranchNames.push_back( labels.at(i) );
      fBranchEntries[labels.at(i)] = {};
    }

    if ( saveRunSubEvt ) {
      assert( fBranchEntries.find("Run/i") == fBranchEntries.end() &&
              fBranchEntries.find("Subrun/i") == fBranchEntries.end() &&
              fBranchEntries.find("Evt/i") == fBranchEntries.end() );

      fOrderedBranchNames.push_back( "Run/i" ); fBranchEntries["Run/i"] = {};
      fOrderedBranchNames.push_back( "Subrun/i" ); fBranchEntries["Subrun/i"] = {};
      fOrderedBranchNames.push_back( "Evt/i" ); fBranchEntries["Evt/i"] = {};
    }

    loader.AddTree( *this, labels, vars, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  // Constructor for a set of MultiVars
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<MultiVar>& vars, const SpillCut& spillcut,
              const Cut& cut, const SystShifts& shift, const bool saveRunSubEvt )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt)
  {
    assert( labels.size() == vars.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      fOrderedBranchNames.push_back( labels.at(i) );
      fBranchEntries[labels.at(i)] = {};
    }

    if ( saveRunSubEvt ) {
      assert( fBranchEntries.find("Run/i") == fBranchEntries.end() &&
              fBranchEntries.find("Subrun/i") == fBranchEntries.end() &&
              fBranchEntries.find("Evt/i") == fBranchEntries.end() );

      fOrderedBranchNames.push_back( "Run/i" ); fBranchEntries["Run/i"] = {};
      fOrderedBranchNames.push_back( "Subrun/i" ); fBranchEntries["Subrun/i"] = {};
      fOrderedBranchNames.push_back( "Evt/i" ); fBranchEntries["Evt/i"] = {};
    }

    loader.AddTree( *this, labels, vars, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  // Constructor for a set of SpillVars
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<SpillVar>& vars, const SpillCut& spillcut, const bool saveRunSubEvt )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt)
  {
    assert( labels.size() == vars.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      fOrderedBranchNames.push_back( labels.at(i) );
      fBranchEntries[labels.at(i)] = {};
    }

    if ( saveRunSubEvt ) {
      assert( fBranchEntries.find("Run/i") == fBranchEntries.end() &&
              fBranchEntries.find("Subrun/i") == fBranchEntries.end() &&
              fBranchEntries.find("Evt/i") == fBranchEntries.end() );

      fOrderedBranchNames.push_back( "Run/i" ); fBranchEntries["Run/i"] = {};
      fOrderedBranchNames.push_back( "Subrun/i" ); fBranchEntries["Subrun/i"] = {};
      fOrderedBranchNames.push_back( "Evt/i" ); fBranchEntries["Evt/i"] = {};
    }

    loader.AddTree( *this, labels, vars, spillcut );
  }

  //----------------------------------------------------------------------
  // Constructor for a set of SpillMultiVars
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<SpillMultiVar>& vars, const SpillCut& spillcut, const bool saveRunSubEvt )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt)
  {
    assert( labels.size() == vars.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      fOrderedBranchNames.push_back( labels.at(i) );
      fBranchEntries[labels.at(i)] = {};
    }

    if ( saveRunSubEvt ) {
      assert( fBranchEntries.find("Run/i") == fBranchEntries.end() &&
              fBranchEntries.find("Subrun/i") == fBranchEntries.end() &&
              fBranchEntries.find("Evt/i") == fBranchEntries.end() );

      fOrderedBranchNames.push_back( "Run/i" ); fBranchEntries["Run/i"] = {};
      fOrderedBranchNames.push_back( "Subrun/i" ); fBranchEntries["Subrun/i"] = {};
      fOrderedBranchNames.push_back( "Evt/i" ); fBranchEntries["Evt/i"] = {};
    }

    loader.AddTree( *this, labels, vars, spillcut );
  }

  //----------------------------------------------------------------------
  // Add an entry to a branch
  void Tree::UpdateEntries( const std::map<std::string, std::vector<double>> valsMap )
  {
    unsigned int idxBranch = 0;
    unsigned int previousSize=0;
    for ( auto const& [name, vals] : valsMap ) {
      if ( idxBranch>0 ) assert(previousSize == vals.size());
      assert ( fBranchEntries.find(name) != fBranchEntries.end() );

      previousSize = vals.size();
      fBranchEntries.at(name).insert(fBranchEntries.at(name).end(),vals.begin(),vals.end());
      idxBranch+=1;
    }
    fNEntries+=(long long)previousSize;
  }

  //----------------------------------------------------------------------
  // Add an entry to a branch
  void Tree::UpdateExposure( const double pot, const double livetime )
  {
    fPOT+=pot;
    fLivetime+=livetime;
  }

  //----------------------------------------------------------------------
  // Save to a ROOT file directory as a TTree -- similar to Spectrum::SaveTo()
  void Tree::SaveTo( TDirectory* dir ) const
  {
    std::cout << "WRITING A TTree FOR THIS Tree OBJECT WITH:" << std::endl;
    std::cout << "  " << fNEntries << " Entries" << std::endl;
    std::cout << "  For " << fPOT << " POT and " << fLivetime << " Livetime" << std::endl;
    if ( fNEntries>=1 ) {
      std::cout << "  Sample:" << std::endl;
      for ( auto const& [key, value] : fBranchEntries ) {
        std::cout << "    " << key << " has first value " << value.at(0) << std::endl;
      }
    }
    else {
      std::cout << "    Double checking --> " << fOrderedBranchNames.at(0)
                << " has size " <<  fBranchEntries.at( fOrderedBranchNames.at(0) ).size() << std::endl;
    }

    // Check (and assert) that the branches all have fNEntries
    for ( auto const& [branch, values] : fBranchEntries ){
      assert( (long long)values.size() == fNEntries );
    }

    TDirectory *tmp = gDirectory;
    dir->cd();

    TH1D thePOT("POT","POT",1,0,1);
    thePOT.SetBinContent(1,fPOT);
    thePOT.Write();

    TH1D theLivetime("Livetime","Livetime",1,0,1);
    theLivetime.SetBinContent(1,fLivetime);
    theLivetime.Write();

    TTree theTree( fTreeName.c_str(), fTreeName.c_str() );

    const int NBranches = fOrderedBranchNames.size();

    bool treatAsInt[ NBranches ];
    double entryValsDouble[ NBranches ];
    long long entryValsInt[ NBranches ];

    for ( unsigned int idxBranch=0; idxBranch<fOrderedBranchNames.size(); ++idxBranch ) {
      if ( fOrderedBranchNames.at(idxBranch).find("/i")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/i")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/I")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/I")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/")!=std::string::npos ) {
        std::cout << "WARNING!! A '/' was found in the variable name, possibly by mistake? Will treat this branch as a double..." << std::endl;
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
      }
      else {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
      }
    }

    // Loop over entries
    for ( unsigned int idxEntry=0; idxEntry < fNEntries; ++idxEntry ) {
      // Fill up the vals
      for ( unsigned int idxBranch=0; idxBranch < fOrderedBranchNames.size(); ++idxBranch ) {
        if ( !treatAsInt[idxBranch] ) entryValsDouble[idxBranch] = fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry);
        else                          entryValsInt[idxBranch] = lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
      }
      theTree.Fill();
    }

    theTree.Write();

    tmp->cd();
  }

}