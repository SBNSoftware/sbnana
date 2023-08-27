#include "sbnana/CAFAna/Core/Tree.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TSpline.h"

#include <cassert>
#include <cmath>

namespace ana
{

  //----------------------------------------------------------------------
  // Constructor for a set of Vars, typical usage for a Selected Tree
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<Var>& vars, const SpillCut& spillcut,
              const Cut& cut, const SystShifts& shift, const bool saveRunSubEvt, const bool saveSliceNum )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(saveSliceNum)
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
    if ( saveSliceNum ) {
      assert( fBranchEntries.find("Slice/i") == fBranchEntries.end() );
      fOrderedBranchNames.push_back( "Slice/i" ); fBranchEntries["Slice/i"] = {};
    }

    loader.AddTree( *this, labels, vars, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  // Constructor for a set of MultiVars
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<MultiVar>& vars, const SpillCut& spillcut,
              const Cut& cut, const SystShifts& shift, const bool saveRunSubEvt, const bool saveSliceNum )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(saveSliceNum)
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
    if ( saveSliceNum ) {
      assert( fBranchEntries.find("Slice/i") == fBranchEntries.end() );
      fOrderedBranchNames.push_back( "Slice/i" ); fBranchEntries["Slice/i"] = {};
    }

    loader.AddTree( *this, labels, vars, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  // Constructor for a set of SpillVars
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<SpillVar>& vars, const SpillCut& spillcut, const bool saveRunSubEvt )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(false)
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
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(false)
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
  // Update branch exposure
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
    bool treatAsLong[ NBranches ];

    double entryValsDouble[ NBranches ];
    int entryValsInt[ NBranches ];
    long long entryValsLong[ NBranches ];

    for ( unsigned int idxBranch=0; idxBranch<fOrderedBranchNames.size(); ++idxBranch ) {
      if ( fOrderedBranchNames.at(idxBranch).find("/i")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/i")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/I")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/I")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/l")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/l")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/L")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/L")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/")!=std::string::npos ) {
        std::cout << "WARNING!! A '/' was found in the variable name, possibly by mistake? Will treat this branch as a double..." << std::endl;
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
      else {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
    }

    // Loop over entries
    for ( unsigned int idxEntry=0; idxEntry < fNEntries; ++idxEntry ) {
      // Fill up the vals
      for ( unsigned int idxBranch=0; idxBranch < fOrderedBranchNames.size(); ++idxBranch ) {
        if ( !treatAsInt[idxBranch] && !treatAsLong[idxBranch] )     entryValsDouble[idxBranch] = fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry);
        else if ( treatAsInt[idxBranch] && !treatAsLong[idxBranch] ) entryValsInt[idxBranch] = (int)lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else if ( treatAsLong[idxBranch] && !treatAsInt[idxBranch] ) entryValsLong[idxBranch] = lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else {
          if ( idxEntry==0 ) std::cout << "ERROR!! Branch " << fOrderedBranchNames.at(idxBranch) << " wants to fill as int and long..." << std::endl;
        }
      }
      theTree.Fill();
    }

    theTree.Write();

    tmp->cd();
  }

  //----------------------------------------------------------------------
  WeightsTree::WeightsTree( const std::string name, const std::vector<std::string>& labels,
                            const unsigned int nSigma, const bool saveRunSubEvt, const bool saveSliceNum, const unsigned int nWeightsExpected )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(saveSliceNum), fNSigma(nSigma), fNWeightsExpected(nWeightsExpected)
  {
    for ( unsigned int i=0; i<labels.size(); ++i ) {
      fOrderedBranchWeightNames.push_back( labels.at(i) );
      fBranchWeightEntries[labels.at(i)] = {};
    }

    if ( saveRunSubEvt ) {
      assert( fBranchEntries.find("Run/i") == fBranchEntries.end() &&
              fBranchEntries.find("Subrun/i") == fBranchEntries.end() &&
              fBranchEntries.find("Evt/i") == fBranchEntries.end() );

      fOrderedBranchNames.push_back( "Run/i" ); fBranchEntries["Run/i"] = {};
      fOrderedBranchNames.push_back( "Subrun/i" ); fBranchEntries["Subrun/i"] = {};
      fOrderedBranchNames.push_back( "Evt/i" ); fBranchEntries["Evt/i"] = {};
    }
    if ( saveSliceNum ) {
      assert( fBranchEntries.find("Slice/i") == fBranchEntries.end() );
      fOrderedBranchNames.push_back( "Slice/i" ); fBranchEntries["Slice/i"] = {};
    }
  }

  //----------------------------------------------------------------------
  void WeightsTree::MergeTree( const Tree& inTree )
  {
    assert( this->fSaveRunSubEvt && this->fSaveSliceNum && inTree.fSaveRunSubEvt && inTree.fSaveSliceNum );

    // This requires the branches to agree and be in the same order...
    // Someone could write a more complicated version but right now this is the foreseen use.
    assert( this->fNEntries == inTree.fNEntries ); // must have the same N Entries
    for ( unsigned int idx = 0; idx < this->fNEntries; ++idx )  {
      // Entries must be the same
      assert( this->fBranchEntries.at("Run/i").at(idx) == inTree.fBranchEntries.at("Run/i").at(idx) &&
              this->fBranchEntries.at("Subrun/i").at(idx) == inTree.fBranchEntries.at("Subrun/i").at(idx) &&
              this->fBranchEntries.at("Evt/i").at(idx) == inTree.fBranchEntries.at("Evt/i").at(idx) &&
              this->fBranchEntries.at("Slice/i").at(idx) == inTree.fBranchEntries.at("Slice/i").at(idx) );
    }

    std::vector<std::string> branchesToFill;
    for ( auto const& branch : inTree.fOrderedBranchNames ) {
      if ( this->fBranchEntries.find( branch ) == this->fBranchEntries.end() ) {
        this->fOrderedBranchNames.push_back( branch );
        this->fBranchEntries[ branch ] = {};
        branchesToFill.push_back( branch );
      }
    }

    // Now do the merge
    for ( unsigned int idx = 0; idx < this->fNEntries; ++idx )  {
      for ( auto const& branch : branchesToFill ) {
        this->fBranchEntries.at(branch).push_back( inTree.fBranchEntries.at(branch).at(idx) );
      }
    }

  }

  //----------------------------------------------------------------------
  void WeightsTree::UpdateEntries( const std::map<std::string, std::vector<double>> valsMap, const std::map<std::string, std::vector<double>> weightMap )
  {
    // First the basic entries (i.e. Run, Subrun, Event or nothing in this case...)
    unsigned int idxBranch = 0;
    unsigned int previousSize=0;
    for ( auto const& [name, vals] : valsMap ) {
      if ( idxBranch>0 ) assert(previousSize == vals.size());
      assert ( fBranchEntries.find(name) != fBranchEntries.end() );
      assert ( vals.size() == 1 ); // this is set up in a way where this cannot be more than one

      previousSize = vals.size();
      fBranchEntries.at(name).insert(fBranchEntries.at(name).end(),vals.begin(),vals.end());
      idxBranch+=1;
    }

    // Now the fNSigma number of weights, stored as a vector per record.
    for ( auto const& [name, weights] : weightMap ) {
      assert ( fBranchWeightEntries.find(name) != fBranchWeightEntries.end() );
      assert ( weights.size() == fNWeightsExpected );

      fBranchWeightEntries.at(name).push_back(weights);
    }

    fNEntries+=1; // we only send 1 record to this at a time in WeightsTree
  }

  //----------------------------------------------------------------------
  void WeightsTree::UpdateExposure( const double pot, const double livetime )
  {
    fPOT+=pot;
    fLivetime+=livetime;
  }

  //----------------------------------------------------------------------
  NSigmasTree::NSigmasTree( const std::string name, const std::vector<std::string>& labels,
                            SpectrumLoaderBase& loader,
                            const std::vector<const ISyst*>& systsToStore, const SpillCut& spillcut,
                            const Cut& cut, const SystShifts& shift, const unsigned int nSigma, const bool saveRunSubEvt, const bool saveSliceNum )
    : WeightsTree(name,labels,nSigma,saveRunSubEvt,saveSliceNum,2*nSigma+1)
  {
    assert( labels.size() == systsToStore.size() );

    loader.AddNSigmasTree( *this, labels, systsToStore, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  void NSigmasTree::SaveToSplines( TDirectory* dir ) const
  {
    std::cout << "WRITING A TTree FOR THIS Tree OBJECT WITH:" << std::endl;
    std::cout << "  " << fNEntries << " Entries" << std::endl;
    std::cout << "  For " << fPOT << " POT and " << fLivetime << " Livetime" << std::endl;
    std::cout << "  Containing " << fOrderedBranchWeightNames.size() << " splines per entry..." << std::endl;

    // Check (and assert) that the branches all have fNEntries
    for ( auto const& [branch, values] : fBranchEntries ){
      assert( (long long)values.size() == fNEntries );
    }
    for ( auto const& [branch, weightVecs] : fBranchWeightEntries ){
      assert( (long long)weightVecs.size() == fNEntries );
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
    bool treatAsLong[ NBranches ];

    double entryValsDouble[ NBranches ];
    int entryValsInt[ NBranches ];
    long long entryValsLong[ NBranches ];

    for ( unsigned int idxBranch=0; idxBranch<fOrderedBranchNames.size(); ++idxBranch ) {
      if ( fOrderedBranchNames.at(idxBranch).find("/i")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/i")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/I")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/I")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/l")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/l")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/L")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/L")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/")!=std::string::npos ) {
        std::cout << "WARNING!! A '/' was found in the variable name, possibly by mistake? Will treat this branch as a double..." << std::endl;
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
      else {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
    }

    const int NBranchesWeights = fOrderedBranchWeightNames.size();
    const int NSigmas = 2*fNSigma+1;

    double sigmasArr[NSigmas];
    for ( unsigned int idxSigma=0; idxSigma<2*fNSigma+1; ++idxSigma ) {
      sigmasArr[idxSigma] = double((-1*int(fNSigma))+int(idxSigma));
    }

    TSpline3 *splinesArr[NBranchesWeights];
    for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
      splinesArr[idxBranchWeight] = nullptr;
      theTree.Branch( fOrderedBranchWeightNames.at(idxBranchWeight).c_str(), &splinesArr[idxBranchWeight] );
    }

    // Loop over entries
    for ( unsigned int idxEntry=0; idxEntry < fNEntries; ++idxEntry ) {
      // Fill up the vals for the standard value branches
      for ( unsigned int idxBranch=0; idxBranch < fOrderedBranchNames.size(); ++idxBranch ) {
        if ( !treatAsInt[idxBranch] && !treatAsLong[idxBranch] )     entryValsDouble[idxBranch] = fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry);
        else if ( treatAsInt[idxBranch] && !treatAsLong[idxBranch] ) entryValsInt[idxBranch] = (int)lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else if ( treatAsLong[idxBranch] && !treatAsInt[idxBranch] ) entryValsLong[idxBranch] = lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else {
          if ( idxEntry==0 ) std::cout << "ERROR!! Branch " << fOrderedBranchNames.at(idxBranch) << " wants to fill as int and long..." << std::endl;
        }
      }
      // Make the splines
      for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
        double weightsArr[NSigmas];
        unsigned int idxVal = 0;
        for ( auto const& val : fBranchWeightEntries.at( fOrderedBranchWeightNames.at(idxBranchWeight) ).at( idxEntry ) ) {
          weightsArr[ idxVal ] = val;
          idxVal+=1;
        }
        TGraph *graph = new TGraph(NSigmas,sigmasArr,weightsArr);
        splinesArr[ idxBranchWeight ] = new TSpline3( TString::Format("%s_%i",fOrderedBranchWeightNames.at(idxBranchWeight).c_str(),idxEntry),
                                                      graph );
      }

      theTree.Fill();
    }

    theTree.Write();

    tmp->cd();
  }

  //----------------------------------------------------------------------
  void NSigmasTree::SaveToGraphs( TDirectory* dir ) const
  {
    std::cout << "WRITING A TTree FOR THIS Tree OBJECT WITH:" << std::endl;
    std::cout << "  " << fNEntries << " Entries" << std::endl;
    std::cout << "  For " << fPOT << " POT and " << fLivetime << " Livetime" << std::endl;
    std::cout << "  Containing " << fOrderedBranchWeightNames.size() << " graphs per entry..." << std::endl;

    // Check (and assert) that the branches all have fNEntries
    for ( auto const& [branch, values] : fBranchEntries ){
      assert( (long long)values.size() == fNEntries );
    }
    for ( auto const& [branch, weightVecs] : fBranchWeightEntries ){
      assert( (long long)weightVecs.size() == fNEntries );
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
    bool treatAsLong[ NBranches ];

    double entryValsDouble[ NBranches ];
    int entryValsInt[ NBranches ];
    long long entryValsLong[ NBranches ];

    for ( unsigned int idxBranch=0; idxBranch<fOrderedBranchNames.size(); ++idxBranch ) {
      if ( fOrderedBranchNames.at(idxBranch).find("/i")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/i")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/I")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/I")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/l")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/l")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/L")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/L")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/")!=std::string::npos ) {
        std::cout << "WARNING!! A '/' was found in the variable name, possibly by mistake? Will treat this branch as a double..." << std::endl;
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
      else {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
    }

    const int NBranchesWeights = fOrderedBranchWeightNames.size();
    const int NSigmas = 2*fNSigma+1;

    double sigmasArr[NSigmas];
    for ( unsigned int idxSigma=0; idxSigma<2*fNSigma+1; ++idxSigma ) {
      sigmasArr[idxSigma] = double((-1*int(fNSigma))+int(idxSigma));
    }

    TGraph *graphsArr[NBranchesWeights];
    for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
      graphsArr[idxBranchWeight] = nullptr;
      theTree.Branch( fOrderedBranchWeightNames.at(idxBranchWeight).c_str(), &graphsArr[idxBranchWeight] );
    }

    // Loop over entries
    for ( unsigned int idxEntry=0; idxEntry < fNEntries; ++idxEntry ) {
      // Fill up the vals for the standard value branches
      for ( unsigned int idxBranch=0; idxBranch < fOrderedBranchNames.size(); ++idxBranch ) {
        if ( !treatAsInt[idxBranch] && !treatAsLong[idxBranch] )     entryValsDouble[idxBranch] = fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry);
        else if ( treatAsInt[idxBranch] && !treatAsLong[idxBranch] ) entryValsInt[idxBranch] = (int)lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else if ( treatAsLong[idxBranch] && !treatAsInt[idxBranch] ) entryValsLong[idxBranch] = lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else {
          if ( idxEntry==0 ) std::cout << "ERROR!! Branch " << fOrderedBranchNames.at(idxBranch) << " wants to fill as int and long..." << std::endl;
        }
      }
      // Make the graphs
      for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
        double weightsArr[NSigmas];
        unsigned int idxVal = 0;
        for ( auto const& val : fBranchWeightEntries.at( fOrderedBranchWeightNames.at(idxBranchWeight) ).at( idxEntry ) ) {
          weightsArr[ idxVal ] = val;
          idxVal+=1;
        }
        graphsArr[ idxBranchWeight ] = new TGraph(NSigmas,sigmasArr,weightsArr);
      }

      theTree.Fill();
    }

    theTree.Write();

    tmp->cd();
  }

  //----------------------------------------------------------------------
  void NSigmasTree::SaveTo( TDirectory* dir ) const
  {
    std::cout << "WRITING A TTree FOR THIS Tree OBJECT WITH:" << std::endl;
    std::cout << "  " << fNEntries << " Entries" << std::endl;
    std::cout << "  For " << fPOT << " POT and " << fLivetime << " Livetime" << std::endl;
    std::cout << "  Containing " << fOrderedBranchWeightNames.size() << " splines per entry..." << std::endl;

    // Check (and assert) that the branches all have fNEntries
    for ( auto const& [branch, values] : fBranchEntries ){
      assert( (long long)values.size() == fNEntries );
    }
    for ( auto const& [branch, weightVecs] : fBranchWeightEntries ){
      assert( (long long)weightVecs.size() == fNEntries );
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
    bool treatAsLong[ NBranches ];

    double entryValsDouble[ NBranches ];
    int entryValsInt[ NBranches ];
    long long entryValsLong[ NBranches ];

    for ( unsigned int idxBranch=0; idxBranch<fOrderedBranchNames.size(); ++idxBranch ) {
      if ( fOrderedBranchNames.at(idxBranch).find("/i")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/i")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/I")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/I")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/l")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/l")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/L")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/L")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/")!=std::string::npos ) {
        std::cout << "WARNING!! A '/' was found in the variable name, possibly by mistake? Will treat this branch as a double..." << std::endl;
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
      else {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
    }

    const unsigned int NSigmas = 2*fNSigma+1;

    std::map< std::string, std::vector<double> > weights;

    // set up map
    for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
      weights[ fOrderedBranchWeightNames.at(idxBranchWeight) ] = {};
    }
    // set up branches
    for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
      theTree.Branch( fOrderedBranchWeightNames.at(idxBranchWeight).c_str(), &weights[fOrderedBranchWeightNames.at(idxBranchWeight)] );
    }

    // Loop over entries
    for ( unsigned int idxEntry=0; idxEntry < fNEntries; ++idxEntry ) {
      // Fill up the vals for the standard value branches
      for ( unsigned int idxBranch=0; idxBranch < fOrderedBranchNames.size(); ++idxBranch ) {
        if ( !treatAsInt[idxBranch] && !treatAsLong[idxBranch] )     entryValsDouble[idxBranch] = fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry);
        else if ( treatAsInt[idxBranch] && !treatAsLong[idxBranch] ) entryValsInt[idxBranch] = (int)lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else if ( treatAsLong[idxBranch] && !treatAsInt[idxBranch] ) entryValsLong[idxBranch] = lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else {
          if ( idxEntry==0 ) std::cout << "ERROR!! Branch " << fOrderedBranchNames.at(idxBranch) << " wants to fill as int and long..." << std::endl;
        }
      }
      // Save the weights
      for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
        for ( unsigned int idxWt = 0; idxWt < NSigmas; ++idxWt ) {
          weights[ fOrderedBranchWeightNames.at(idxBranchWeight) ].push_back(fBranchWeightEntries.at( fOrderedBranchWeightNames.at(idxBranchWeight) ).at(idxEntry).at(idxWt));
        }
      }

      theTree.Fill();

      // clear eights vec:
      for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
        weights[ fOrderedBranchWeightNames.at(idxBranchWeight) ].clear();
      }
    }

    theTree.Write();

    tmp->cd();
  }

  //----------------------------------------------------------------------
  NUniversesTree::NUniversesTree( const std::string name, const std::vector<std::string>& labels,
                            SpectrumLoaderBase& loader,
                            const std::vector<std::vector<Var>>& univsKnobs, const unsigned int nUniverses,
                            const SpillCut& spillcut,
                            const Cut& cut, const SystShifts& shift, const bool saveRunSubEvt, const bool saveSliceNum )
    : WeightsTree(name,labels,nUniverses,saveRunSubEvt,saveSliceNum,nUniverses)
  {
    assert( labels.size() == univsKnobs.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      assert( univsKnobs.at(i).size() == fNSigma );
    }

    loader.AddNUniversesTree( *this, labels, univsKnobs, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  void NUniversesTree::SaveTo( TDirectory* dir ) const
  {
    std::cout << "WRITING A TTree FOR THIS Tree OBJECT WITH:" << std::endl;
    std::cout << "  " << fNEntries << " Entries" << std::endl;
    std::cout << "  For " << fPOT << " POT and " << fLivetime << " Livetime" << std::endl;
    std::cout << "  Containing " << fOrderedBranchWeightNames.size() << " splines per entry..." << std::endl;

    // Check (and assert) that the branches all have fNEntries
    for ( auto const& [branch, values] : fBranchEntries ){
      assert( (long long)values.size() == fNEntries );
    }
    for ( auto const& [branch, weightVecs] : fBranchWeightEntries ){
      assert( (long long)weightVecs.size() == fNEntries );
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
    bool treatAsLong[ NBranches ];

    double entryValsDouble[ NBranches ];
    int entryValsInt[ NBranches ];
    long long entryValsLong[ NBranches ];

    for ( unsigned int idxBranch=0; idxBranch<fOrderedBranchNames.size(); ++idxBranch ) {
      if ( fOrderedBranchNames.at(idxBranch).find("/i")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/i")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/I")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/I")).c_str(),
                        &entryValsInt[idxBranch] );
        treatAsInt[idxBranch] = true;
        treatAsLong[idxBranch] = false;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/l")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/l")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/L")!=std::string::npos ) {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).substr(0, fOrderedBranchNames.at(idxBranch).find("/L")).c_str(),
                        &entryValsLong[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = true;
      }
      else if ( fOrderedBranchNames.at(idxBranch).find("/")!=std::string::npos ) {
        std::cout << "WARNING!! A '/' was found in the variable name, possibly by mistake? Will treat this branch as a double..." << std::endl;
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
      else {
        theTree.Branch( fOrderedBranchNames.at(idxBranch).c_str(), &entryValsDouble[idxBranch] );
        treatAsInt[idxBranch] = false;
        treatAsLong[idxBranch] = false;
      }
    }

    std::map< std::string, std::vector<double> > weights;

    // set up map
    for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
      weights[ fOrderedBranchWeightNames.at(idxBranchWeight) ] = {};
    }
    // set up branches
    for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
      theTree.Branch( fOrderedBranchWeightNames.at(idxBranchWeight).c_str(), &weights[fOrderedBranchWeightNames.at(idxBranchWeight)] );
    }

    // Loop over entries
    for ( unsigned int idxEntry=0; idxEntry < fNEntries; ++idxEntry ) {
      // Fill up the vals for the standard value branches
      for ( unsigned int idxBranch=0; idxBranch < fOrderedBranchNames.size(); ++idxBranch ) {
        if ( !treatAsInt[idxBranch] && !treatAsLong[idxBranch] )     entryValsDouble[idxBranch] = fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry);
        else if ( treatAsInt[idxBranch] && !treatAsLong[idxBranch] ) entryValsInt[idxBranch] = (int)lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else if ( treatAsLong[idxBranch] && !treatAsInt[idxBranch] ) entryValsLong[idxBranch] = lround(fBranchEntries.at( fOrderedBranchNames.at(idxBranch) ).at(idxEntry));
        else {
          if ( idxEntry==0 ) std::cout << "ERROR!! Branch " << fOrderedBranchNames.at(idxBranch) << " wants to fill as int and long..." << std::endl;
        }
      }
      // Save the weights
      for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
        for ( unsigned int idxWt = 0; idxWt < fNSigma; ++idxWt ) {
          weights[ fOrderedBranchWeightNames.at(idxBranchWeight) ].push_back(fBranchWeightEntries.at( fOrderedBranchWeightNames.at(idxBranchWeight) ).at(idxEntry).at(idxWt));
        }
      }

      theTree.Fill();

      // clear eights vec:
      for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
        weights[ fOrderedBranchWeightNames.at(idxBranchWeight) ].clear();
      }
    }

    theTree.Write();

    tmp->cd();
  }

}
