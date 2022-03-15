#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "SBNAna/Vars/Vars.h"
#include "SBNAna/Vars/NumuVars.h"
#include "SBNAna/Cuts/NumuCuts.h"
#include "SBNAna/Cuts/TruthCuts.h"

using namespace ana;


// Useful for counting the number of events that pass a cut
const Var kCounting([](const caf::SRSliceProxy *slc)
		    {
		      return 1;
		    });

// Costh of the numu primary track (See Vars/NumuVars.cxx)
const Var kPrimTrkCosth([](const caf::SRSliceProxy *slc)
			{
			  int muIdx = kPrimMuonIdx(slc);
			  if( muIdx < 0 ) return -5.;

			  return (double)slc->reco.trk[muIdx].costh;
			});

// Costh of the numu primary track (See Vars/NumuVars.cxx)
const Var kPrimTrkCRTdist([](const caf::SRSliceProxy *slc)
			{
			  int muIdx = kPrimMuonIdx(slc);
			  if( muIdx < 0 ) return -5.;
			  
			  return (double)slc->reco.trk[muIdx].crthit.distance;
			});


// These are examples of useful structs to
// use when making a bunch of Spectra
struct PlotDef
{
  std::string suffix = "";
  std::string label = "";
  Binning bins = Binning::Simple(3,0,3);
  Var var = kCounting;
};

// In this example, we are making the following Spectra
std::vector<PlotDef> plots = 
  {{"count",   "",          Binning::Simple(3,0,3),       kCounting},
   {"mucosth", "cos#theta", Binning::Simple(50,-1,1),     kPrimTrkCosth},
   {"vtxx",    "X (cm)",    Binning::Simple(50,-200,200), kSlcVtxX},
   {"vtxy",    "Y (cm)",    Binning::Simple(50,-200,200), kSlcVtxY},
   {"vtxz",    "Z (cm)",    Binning::Simple(50,0,500),    kSlcVtxZ},
   {"mulen", "Primary Track Length", Binning::Simple(50,0,200),     kPrimTrkLen},
   {"mucrtdist", "Primary Track Distance to CRT hit", Binning::Simple(50,0,100),     kPrimTrkCRTdist},
  };

// Selection Struc
struct SelDef
{
  std::string suffix = "";
  std::string label = "";
  Cut cut = kNoCut;
  int color = kBlack;
};

// In this example, our signal is Numu cc
// See Cuts/TruthCuts to check out these cuts
const Cut kSig = kIsNumu && !kIsNC;

// We are making the Spectra defined above 
// for 3 different selections. 
std::vector<SelDef> sels = 
  {{"nocut", "All Slices",  kNoCut, kBlack},
   {"sig",   "True NumuCC", kSig,   kRed+1},
   {"bkg",   "Not NumuCC",  !kSig,  kAzure+2},
  };

// If you wanted to  add a new cut, define an
// expression for kYourCut (see examples in SBNAna/Cuts)
// and then add the following lines to the sels vector: 
/* {"all_yourcut", "All Slices",  kYourCut,            kBlack}, */
/* {"sig_yourcut", "True NumuCC", (kSig && kYourCut),  kRed+1}, */
/* {"bkg_yourcut", "Not NumuCC",  (!kSig && kYourCut), kAzure+2}, */
