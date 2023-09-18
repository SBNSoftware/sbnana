#include "sbnana/CAFAna/Core/Tree.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TClonesArray.h"

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
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(saveSliceNum), fSaveTruthCutType(false), SignalSelection(kNoCut)
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
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(saveSliceNum), fSaveTruthCutType(false), SignalSelection(kNoCut)
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
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(false), fSaveTruthCutType(false), SignalSelection(kNoCut)
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
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(false), fSaveTruthCutType(false), SignalSelection(kNoCut)
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
  // Constructor for a set of TruthVars
  Tree::Tree( const std::string name, const std::vector<std::string>& labels,
              SpectrumLoaderBase& loader,
              const std::vector<TruthVar>& vars, const SpillCut& spillcut,
              const TruthCut& truthcut,
              const Cut& SignalSelection,
              const SystShifts& shift,
              const bool saveRunSubEvt)
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(false), fSaveTruthCutType(true), SignalSelection(SignalSelection)
  {
    assert( labels.size() == vars.size() );

    fOrderedBranchNames.push_back( "CutType/i" ); fBranchEntries["CutType/i"] = {};
    fOrderedBranchNames.push_back( "SpillCutType/i" ); fBranchEntries["SpillCutType/i"] = {};

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

    loader.AddTree( *this, labels, vars, spillcut, truthcut, shift );
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
  CovarianceMatrixTree::CovarianceMatrixTree( const std::string name, const std::vector<std::string>& labels, const std::string weightLabel,
                                              const std::vector< std::vector<std::pair<std::string,std::pair<double,double>>> >& binning,
                                              SpectrumLoaderBase& loader, const std::vector<SpillMultiVar>& vars, const unsigned int nUniverses,
                                              const SpillCut& spillcut )
    : Tree(name,labels,loader,vars,spillcut,false),fWeightLabel(weightLabel),fNUniverses(nUniverses)
  {
    assert( labels.size() == vars.size() );
    fWeights = {};

    // Reconstruct a binning map from the passed vector:
    unsigned int idxBin = 1;
    for ( auto const& binCuts : binning ) {
      fBinning[idxBin] = binCuts;
      idxBin+=1;
    }

    // Check that the cuts make sense:
    // Logic: any cut that appears in the first bin must appear in the rest, and
    //        nothing can appear in a latter bin that isn't in the first
    // TODO: we DO NOT check that no two bins have the same cuts...
    // TODO: also we DO NOT check, but assume later, that if it's a one-value bin
    //       in the first case then it's always this way...
    std::set<std::string> tmpCuts;

    for ( auto const& [bin, cuts] : fBinning ) {
      if ( tmpCuts.size()==0 ) {
        for ( auto const& cut : cuts ) {
          assert(tmpCuts.find(cut.first)==tmpCuts.end());
          tmpCuts.emplace( cut.first );
        }
      }
      else {
        assert(cuts.size()==tmpCuts.size());
        for ( auto const& cut : cuts ) {
          assert(tmpCuts.find(cut.first)!=tmpCuts.end());
        }
      }
    }
    assert(tmpCuts.size()!=0);
  }

  //----------------------------------------------------------------------
  void CovarianceMatrixTree::UpdateEntries( const std::map<std::string, std::vector<double>> valsMap )
  {
    unsigned int idxBranch = 0;
    unsigned int previousSize=0;
    for ( auto const& [name, vals] : valsMap ) {
      if ( idxBranch>0 && name.compare(fWeightLabel)!=0 ) assert(previousSize == vals.size());
      else if ( idxBranch>0 && name.compare(fWeightLabel)==0 ) assert(previousSize == vals.size()/fNUniverses);

      if ( name.compare(fWeightLabel)!=0 ) {
        assert ( fBranchEntries.find(name) != fBranchEntries.end() );
        previousSize = vals.size();
        fBranchEntries.at(name).insert(fBranchEntries.at(name).end(),vals.begin(),vals.end());
        idxBranch+=1;
      }
      else {
        previousSize = vals.size()/fNUniverses;
        std::vector<double> thisWts;        
        for ( unsigned int idxWt=0; idxWt < vals.size(); ++idxWt ) {
          thisWts.push_back( vals.at(idxWt) );
          if ( thisWts.size() == fNUniverses ) {
            fWeights.push_back(thisWts);
            thisWts.clear();
          }
        }
        // Still fill the main branch
        fBranchEntries.at(name).insert(fBranchEntries.at(name).end(),vals.begin(),vals.end());
      }
    }

    fNEntries+=(long long)previousSize;
  }

  //----------------------------------------------------------------------
  void CovarianceMatrixTree::SaveTo( TDirectory* dir ) const
  {
    assert(fOrderedBranchNames.size() >= 2); // need at least a weight label and a key for the binning...

    // Pick a random key and double check the sizes make sense
    std::string key = fOrderedBranchNames.at(0);
    if ( key.compare(fWeightLabel)==0 ) key = fOrderedBranchNames.at(1);
    assert( fNEntries==(long long)fWeights.size() && fBranchEntries.at(key).size()==fWeights.size() );

    std::cout << "SAVING a covariance matrix for the variable labeled " << fWeightLabel << " with binning:" << std::endl;
    PrintBinning();

    TDirectory *tmp = gDirectory;
    dir->cd();

    TH1D thePOT("POT","POT",1,0,1);
    thePOT.SetBinContent(1,fPOT);
    thePOT.Write();

    TH1D theLivetime("Livetime","Livetime",1,0,1);
    theLivetime.SetBinContent(1,fLivetime);
    theLivetime.Write();

    // Make the covariance matrix
    const unsigned int nBins = fBinning.size();

    // construct the bins to sum weights
    std::vector<std::vector<double>> binCts;
    for ( unsigned int idxU = 0; idxU < fNUniverses; ++idxU ) {
      std::vector<double> binCtsU;
      for ( unsigned int idxB = 0; idxB < nBins; ++idxB ) {
        binCtsU.push_back(0.);
      }
      binCts.push_back( binCtsU );
    }

    // fill them 
    for ( unsigned int idxEntry = 0; idxEntry < fNEntries; ++idxEntry ) {
      // determine bin:
      int theBin = 0;
      for ( auto const& [bin, cuts] : fBinning ) {
        bool failedCut = false;
        for ( auto const& cut : cuts ) {
          if ( failedCut ) continue;
          assert( fBranchEntries.find(cut.first)!=fBranchEntries.end() );
          const double loVal = cut.second.first;
          const double hiVal = cut.second.second;
          // Treat as long ints if needed
          if ( cut.first.find("/i")!=std::string::npos ||
               cut.first.find("/I")!=std::string::npos ) {
            const int loValI = std::lround(loVal);
            const int hiValI = std::lround(hiVal);
            const int thisVal = std::lround( fBranchEntries.at(cut.first).at(idxEntry) );
            if ( loValI==hiValI ) {
              if ( thisVal!=loValI ) {
                failedCut = true;
              }
            }
            else if ( thisVal <= loValI || thisVal > hiValI ){
              failedCut = true;
            }
          }
          else if ( cut.first.find("/l")!=std::string::npos ||
                    cut.first.find("/L")!=std::string::npos ) {
            const long long loValI = std::lround(loVal);
            const long long hiValI = std::lround(hiVal);
            const long long thisVal = std::lround( fBranchEntries.at(cut.first).at(idxEntry) );
            if ( loValI==hiValI ) {
              if ( thisVal!=loValI ) {
                failedCut = true;
              }
            }
            else if ( thisVal <= loValI || thisVal > hiValI ){
              failedCut = true;
            }
          }
          // double check if the lo and hi are the same
          else if ( fabs(hiVal-loVal) < std::numeric_limits<double>::epsilon() ) {
            std::cout << "WARNING!!! You may have been expecting " << cut.first << " to be an integer but it's being treated like a double!" << std::endl;
            const double thisVal = fBranchEntries.at(cut.first).at(idxEntry);
            if ( fabs(thisVal-loVal) > std::numeric_limits<double>::epsilon() ) failedCut = true;
          }
          // else check the range, inclusive on high edge
          else {
            const double thisVal = fBranchEntries.at(cut.first).at(idxEntry);
            if ( thisVal <= loVal || thisVal > hiVal ){
              failedCut = true;
            }
          }
        }
        // if got to the end of cuts and we pass, then this is our bin
        if( !failedCut ){
          theBin = bin;
          break;
        }
      }
      // if under- or over-flow then skip this event
      if ( theBin == 0 ) continue;

      // fill for each universe:
      for ( unsigned int idxU = 0; idxU < fNUniverses; ++idxU ) {
        //if ( idxEntry < 10 ) std::cout << "idxEntry " << idxEntry << ": filling bin " << theBin << " with value " << (fWeights.at(idxEntry).at(idxU)-1.) << " for universe " << idxU << std::endl;
        binCts.at(idxU).at(theBin-1)+=(fWeights.at(idxEntry).at(idxU)); // Was thinking to subtract 1 to make it relative contributions but it's all relative to the mean so I think I don't want to do that...
      }
    }

    // Make the vector of means
    std::vector<double> means;
    for ( unsigned int idxB = 0; idxB < nBins; ++idxB ) {
      double sum = 0.;
      for ( unsigned int idxU = 0; idxU < fNUniverses; ++idxU ){
        sum+=binCts.at(idxU).at(idxB);
      }

      means.push_back( sum / fNUniverses );
    }

    // To make the correlation matrix we also want standard deviations
    std::vector<double> sigmas;
    for ( unsigned int idxB = 0; idxB < nBins; ++idxB ) {
      double sum = 0.;
      for ( unsigned int idxU = 0; idxU < fNUniverses; ++idxU ){
        sum+=std::pow(binCts.at(idxU).at(idxB)-means.at(idxB),2);
      }

      sigmas.push_back( std::sqrt(sum / (fNUniverses-1)) );
    }

    //std::cout << "Means and sigmas per bin:" << std::endl;
    //for ( unsigned int idxB = 0; idxB < nBins; ++idxB ) {
    //  std::cout << "Bin " << idxB+1 << ", mean: " << means.at(idxB) << ", sigma: " << sigmas.at(idxB) << std::endl;
    //}

    // Now make a TH2D to write the covariance (using formula as in Forumula 10 of sbn-doc-27384):
    TH2D *h_cov_mtx = new TH2D( TString::Format("hCovMtx_%s",fWeightLabel.c_str()), ";bins;bins", nBins, 1, nBins+1, nBins, 1, nBins+1 );

    for ( unsigned int xBin=1; xBin<=nBins; ++xBin ) {
      for ( unsigned int yBin=1; yBin<=nBins; ++yBin ) {
        double cov = 0.;
        for ( unsigned int idxU = 0; idxU < fNUniverses; ++idxU ) {
          //if ( xBin < 3 ) std::cout << "xBin " << xBin << ", yBin " << yBin << ": Adding to covariance " << ((binCts.at(idxU).at(xBin-1)-means.at(xBin-1))*(binCts.at(idxU).at(yBin-1)-means.at(yBin-1))) << " for universe " << idxU << std::endl;
          cov+=((binCts.at(idxU).at(xBin-1)-means.at(xBin-1))*(binCts.at(idxU).at(yBin-1)-means.at(yBin-1)));
        }
        cov/=(fNUniverses-1);

        h_cov_mtx->SetBinContent(xBin,yBin,cov);
      }
    }

    h_cov_mtx->Write();

    // and using the documentation in https://root.cern/doc/master/group__Matrix.html
    // Let's also save this as a TMatrixD
    // TODO: do I need to save overflow/underflow?
    TArrayD tarray(nBins*nBins);

    for ( unsigned int xBin=1; xBin<=nBins; ++xBin ) {
      const unsigned int row = xBin-1;
      for ( unsigned int yBin=1; yBin<=nBins; ++yBin ) {
        const unsigned int column = yBin-1;
        const double val = h_cov_mtx->GetBinContent(xBin,yBin);
        tarray[nBins*row+column] = val;
      }
    }

    TMatrixD cov_mtx(nBins,nBins);
    cov_mtx.SetMatrixArray(tarray.GetArray());

    cov_mtx.Write( TString::Format("hCovMtx_%s",fWeightLabel.c_str()) );

    // And correlation matrix, by way of Corr(X,Y) = Cov(X,Y)/Sigma(X)Sigma(Y)
    // See e.g. PDG 2018, Chapter 39 Page 25, or https://en.wikipedia.org/wiki/Correlation#Correlation_matrices
    TH2D *h_corr_mtx = new TH2D( TString::Format("hCorrMtx_%s",fWeightLabel.c_str()), ";bins;bins", nBins, 1, nBins+1, nBins, 1, nBins+1 );

    for ( unsigned int xBin=1; xBin<=nBins; ++xBin ) {
      for ( unsigned int yBin=1; yBin<=nBins; ++yBin ) {
        double corr = h_cov_mtx->GetBinContent(xBin,yBin) / (sigmas[xBin-1]*sigmas[yBin-1]);

        h_corr_mtx->SetBinContent(xBin,yBin,corr);
      }
    }

    h_corr_mtx->Write();

    TArrayD tarray_corr(nBins*nBins);
    for ( unsigned int xBin=1; xBin<=nBins; ++xBin ) {
      const unsigned int row = xBin-1;
      for ( unsigned int yBin=1; yBin<=nBins; ++yBin ) {
        const unsigned int column = yBin-1;
        const double val = h_corr_mtx->GetBinContent(xBin,yBin);
        tarray_corr[nBins*row+column] = val;
      }
    }
    TMatrixD corr_mtx(nBins,nBins);
    corr_mtx.SetMatrixArray(tarray_corr.GetArray());

    corr_mtx.Write( TString::Format("hCorrMtx_%s",fWeightLabel.c_str()) );

    tmp->cd();
  }

  //----------------------------------------------------------------------
  void CovarianceMatrixTree::PrintBinning() const
  {
    std::vector<std::string> cutNames;
    std::vector<bool> treatAsInt;
    std::vector<bool> treatAsLong;
    std::vector<bool> treatOnlyLo;

    for ( auto const& [bin, cuts] : fBinning ) {
      if ( cutNames.size()==0 ) {
        std::string HdrString = "variables: ";
        for ( auto const& cut : cuts ) {
          cutNames.push_back( cut.first );

          if ( cut.first.find("/i")!=std::string::npos ||
               cut.first.find("/I")!=std::string::npos ) {
            treatAsInt.push_back(true);
            treatAsLong.push_back(false);
            int checkLo = std::lround(cut.second.first);
            int checkHi = std::lround(cut.second.second);
            if ( checkLo==checkHi ) treatOnlyLo.push_back(true);
            else treatOnlyLo.push_back(false);

            HdrString+=(cut.first+" ");
          }
          else if ( cut.first.find("/l")!=std::string::npos ||
                    cut.first.find("/L")!=std::string::npos ) {
            treatAsInt.push_back(false);
            treatAsLong.push_back(true);
            long long checkLo = std::lround(cut.second.first);
            long long checkHi = std::lround(cut.second.second);
            if ( checkLo==checkHi ) treatOnlyLo.push_back(true);
            else treatOnlyLo.push_back(false);

            HdrString+=(cut.first+" ");
          }
          else {
            treatAsInt.push_back(false);
            treatAsLong.push_back(false);
            treatOnlyLo.push_back(false);
            HdrString+=(cut.first+" "+cut.first+" ");
          }
        }
        std::cout << HdrString << std::endl;
      }

      std::string BinString;
      for ( unsigned int idxCut = 0; idxCut<cutNames.size(); ++idxCut ) {
        auto const& cut = fBinning.at(bin).at(idxCut);
        const double loVal = cut.second.first;
        const double hiVal = cut.second.second;

        if ( treatAsInt[idxCut] ) { 
          const int loValI = std::lround(loVal);
          const int hiValI = std::lround(hiVal);

          if ( treatOnlyLo[idxCut] ) BinString+=(std::to_string(loValI)+" ");
          else BinString+=(std::to_string(loValI)+" "+std::to_string(hiValI)+" ");
        }
        else if ( treatAsLong[idxCut] ) { 
          const long long loValI = std::lround(loVal);
          const long long hiValI = std::lround(hiVal);

          if ( treatOnlyLo[idxCut] ) BinString+=(std::to_string(loValI)+" ");
          else BinString+=(std::to_string(loValI)+" "+std::to_string(hiValI)+" ");
        }
        else BinString+=(std::to_string(loVal)+" "+std::to_string(hiVal)+" ");
      }
      std::cout << BinString << std::endl;

    }
  }

  //----------------------------------------------------------------------
  WeightsTree::WeightsTree( const std::string name, const std::vector<std::string>& labels,
                            const std::vector<unsigned int>& nWeights, const bool saveRunSubEvt, const bool saveSliceNum )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(saveSliceNum)
  {
    assert( nWeights.size()==labels.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      fOrderedBranchWeightNames.push_back( labels.at(i) );
      fBranchWeightEntries[labels.at(i)] = {};

      fNSigmasLo[labels.at(i)] = 0;
      fNSigmasHi[labels.at(i)] = nWeights.at(i);
      fNWeightsExpected[labels.at(i)] = nWeights.at(i);
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
  WeightsTree::WeightsTree( const std::string name, const std::vector<std::string>& labels,
                            const std::vector<std::pair<int,int>>& nSigma, const bool saveRunSubEvt, const bool saveSliceNum )
    : fTreeName(name), fNEntries(0), fPOT(0), fLivetime(0), fSaveRunSubEvt(saveRunSubEvt), fSaveSliceNum(saveSliceNum)
  {
    assert( nSigma.size()==labels.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      fOrderedBranchWeightNames.push_back( labels.at(i) );
      fBranchWeightEntries[labels.at(i)] = {};

      assert( nSigma.at(i).second > nSigma.at(i).first );

      fNSigmasLo[labels.at(i)] = nSigma.at(i).first;
      fNSigmasHi[labels.at(i)] = nSigma.at(i).second;
      fNWeightsExpected[labels.at(i)] = (unsigned int)((nSigma.at(i).second-nSigma.at(i).first) + 1);
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
    // RunSubEvt must be saved
    assert( this->fSaveRunSubEvt && inTree.fSaveRunSubEvt );
    // Tree with SpilVar or TruthVar does not require SliceNum to be saved.
    // But if either of "This" or "inTree" does save SliceNum, then both should have it saved
    bool SliceNumIsSaved = false;
    if( this->fSaveSliceNum || inTree.fSaveSliceNum ){
      assert( this->fSaveSliceNum && inTree.fSaveSliceNum );
      SliceNumIsSaved = true;
    }

    // This requires the branches to agree and be in the same order...
    // Someone could write a more complicated version but right now this is the foreseen use.
    assert( this->fNEntries == inTree.fNEntries ); // must have the same N Entries
    for ( unsigned int idx = 0; idx < this->fNEntries; ++idx )  {
      // Entries must be the same
      assert( this->fBranchEntries.at("Run/i").at(idx) == inTree.fBranchEntries.at("Run/i").at(idx) &&
              this->fBranchEntries.at("Subrun/i").at(idx) == inTree.fBranchEntries.at("Subrun/i").at(idx) &&
              this->fBranchEntries.at("Evt/i").at(idx) == inTree.fBranchEntries.at("Evt/i").at(idx) );
      if(SliceNumIsSaved){
        assert( this->fBranchEntries.at("Slice/i").at(idx) == inTree.fBranchEntries.at("Slice/i").at(idx) );
      }
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

    // Now the fNWeightsExpected number of weights, stored as a vector per record.
    for ( auto const& [name, weights] : weightMap ) {
      assert ( fBranchWeightEntries.find(name) != fBranchWeightEntries.end() );
      assert ( weights.size() == fNWeightsExpected.at(name) );

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
                            const std::vector<const ISyst*>& systsToStore, const std::vector<std::pair<int,int>>& nSigma,
                            const SpillCut& spillcut,
                            const Cut& cut, const SystShifts& shift, const bool saveRunSubEvt, const bool saveSliceNum )
  : WeightsTree(name,labels,nSigma,saveRunSubEvt,saveSliceNum)
  {
    assert( labels.size() == systsToStore.size() );
    assert( nSigma.size() == labels.size() );

    loader.AddNSigmasTree( *this, labels, systsToStore, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  NSigmasTree::NSigmasTree( const std::string name, const std::vector<std::string>& labels,
                            SpectrumLoaderBase& loader,
                            const std::vector<const ISyst*>& systsToStore, const std::vector<std::pair<int,int>>& nSigma,
                            const TruthCut& truthcut, const SystShifts& shift, const bool saveRunSubEvt)
  : WeightsTree(name,labels,nSigma,saveRunSubEvt,false)
  {
    assert( labels.size() == systsToStore.size() );
    assert( nSigma.size() == labels.size() );

    loader.AddNSigmasTree( *this, labels, systsToStore, truthcut, shift );
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
        const int NSigmas = fNWeightsExpected.at( fOrderedBranchWeightNames.at(idxBranchWeight) );
        double sigmasArr[NSigmas];
        for ( unsigned int idxSigma=0; idxSigma<(unsigned int)NSigmas; ++idxSigma ) {
          sigmasArr[idxSigma] = double(fNSigmasLo.at( fOrderedBranchWeightNames.at(idxBranchWeight) ) + int(idxSigma));
        }

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
        const int NSigmas = fNWeightsExpected.at( fOrderedBranchWeightNames.at(idxBranchWeight) );
        double sigmasArr[NSigmas];
        for ( unsigned int idxSigma=0; idxSigma<(unsigned int)NSigmas; ++idxSigma ) {
          sigmasArr[idxSigma] = double(fNSigmasLo.at( fOrderedBranchWeightNames.at(idxBranchWeight) ) + int(idxSigma));
        }

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
  void NSigmasTree::SaveToTClonesArrays( TDirectory* dir ) const
  {

    std::cout << "WRITING A TTree FOR THIS Tree OBJECT WITH:" << std::endl;
    std::cout << "  " << fNEntries << " Entries" << std::endl;
    std::cout << "  For " << fPOT << " POT and " << fLivetime << " Livetime" << std::endl;
    std::cout << "  Containing " << fOrderedBranchWeightNames.size() << " TClonesArray per entry..." << std::endl;

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

    TClonesArray *arraysArr[NBranchesWeights];
    for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
      arraysArr[idxBranchWeight] = new TClonesArray("TGraph", 1);
      theTree.Branch( fOrderedBranchWeightNames.at(idxBranchWeight).c_str(), &arraysArr[idxBranchWeight], 32000, -1);
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
        const int NSigmas = fNWeightsExpected.at( fOrderedBranchWeightNames.at(idxBranchWeight) );
        double sigmasArr[NSigmas];
        for ( unsigned int idxSigma=0; idxSigma<(unsigned int)NSigmas; ++idxSigma ) {
          sigmasArr[idxSigma] = double(fNSigmasLo.at( fOrderedBranchWeightNames.at(idxBranchWeight) ) + int(idxSigma));
        }

        double weightsArr[NSigmas];
        unsigned int idxVal = 0;
        for ( auto const& val : fBranchWeightEntries.at( fOrderedBranchWeightNames.at(idxBranchWeight) ).at( idxEntry ) ) {
          weightsArr[ idxVal ] = val;
          idxVal+=1;
        }

        new( (*arraysArr[idxBranchWeight])[0]) TGraph(NSigmas,sigmasArr,weightsArr);

      }

      theTree.Fill();

      for ( unsigned int idxBranchWeight=0; idxBranchWeight<fOrderedBranchWeightNames.size(); ++idxBranchWeight ) {
        arraysArr[idxBranchWeight]->Clear();
      }

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
        const int NSigmas = fNWeightsExpected.at( fOrderedBranchWeightNames.at(idxBranchWeight) );
        for ( unsigned int idxWt = 0; idxWt < (unsigned int)NSigmas; ++idxWt ) {
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
                            const std::vector<std::vector<Var>>& univsKnobs, const std::vector<unsigned int>& nUniverses,
                            const SpillCut& spillcut,
                            const Cut& cut, const SystShifts& shift, const bool saveRunSubEvt, const bool saveSliceNum )
  : WeightsTree(name,labels,nUniverses,saveRunSubEvt,saveSliceNum)
  {
    assert( labels.size() == univsKnobs.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      assert( univsKnobs.at(i).size() == nUniverses.at(i) );
    }

    loader.AddNUniversesTree( *this, labels, univsKnobs, spillcut, cut, shift );
  }

  //----------------------------------------------------------------------
  NUniversesTree::NUniversesTree( const std::string name, const std::vector<std::string>& labels,
                                  SpectrumLoaderBase& loader,
                                  const std::vector<std::vector<TruthVar>>& univsKnobs,
                                  const std::vector<unsigned int>& nUniverses,
                                  const TruthCut& truthcut,
                                  const SystShifts& shift, const bool saveRunSubEvt)
  : WeightsTree(name,labels,nUniverses,saveRunSubEvt,false)
  {
    assert( labels.size() == univsKnobs.size() );

    for ( unsigned int i=0; i<labels.size(); ++i ) {
      assert( univsKnobs.at(i).size() == nUniverses.at(i) );
    }

    loader.AddNUniversesTree( *this, labels, univsKnobs, truthcut, shift );
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
        for ( unsigned int idxWt = 0; idxWt < fNWeightsExpected.at( fOrderedBranchWeightNames.at(idxBranchWeight) ); ++idxWt ) {
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
