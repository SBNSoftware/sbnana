#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/SBNAna/Vars/NumuVarsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/NumuCuts.h"
#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"

using namespace ana;


// Useful for counting the number of events that pass a cut
const Var kCounting_scale([](const caf::SRSliceProxy *slc)
		    {
                      double count = 6e20/3.08823e19;
		      return count;
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
   {"count_scale",   "",    Binning::Simple(3,0,3),       kCounting_scale},
   {"mucosth", "cos#theta", Binning::Simple(50,-1,1),     kPrimTrkCosth},
   {"vtxx",    "X (cm)",    Binning::Simple(50,-200,200), kSlcVtxX},
   {"vtxy",    "Y (cm)",    Binning::Simple(50,-200,200), kSlcVtxY},
   {"vtxz",    "Z (cm)",    Binning::Simple(50,0,500),    kSlcVtxZ},
   {"mulen", "Primary Track Length", Binning::Simple(50,0,200),     kPrimTrkLen},
   {"mucrtdist", "Primary Track Distance to CRT hit", Binning::Simple(50,0,100),     kPrimTrkCRTdist},
   {"CRTangle", "CRT Track Angle", Binning::Simple(20,0,1),     kCRTTrkAngle},
   {"CRTHitTime", "CRT Hit Time", Binning::Simple(20,0,100),     kCRTTrkTime},
   {"CRTHitDist", "CRT Hit Distance", Binning::Simple(20,0,100),     kCRTHitDist},
   {"muonprimtrkp", "Primary Track Momentum", Binning::Simple(24,0,2),     kPrimaryMuonTrkP},
   {"muonprimtrklen", "Primary Track Lenght", Binning::Simple(18,40,400),     kPrimaryMuonTrkLen},
   //{"flashmatch", "Flash Score", Binning::Simple(25,0,25),     kSlcHasFlashMatch},
   {"nuscore", "Pandora Nu Score", Binning::Simple(25,0,1),     kNuScore},
  };

// Selection Struc
struct SelDef
{
  std::string suffix = "";
  std::string label = "";
  Cut cut = kNoCut;
  int color = kBlack;
};

//const Cut kIsNuSlice = ( kTruthIndex >= 0.f );
//const Cut kIsNuSlice = kIsNumu && !kIsNC;
const Cut kIsNuSlice([](const caf::SRSliceProxy* slc) {
      return (slc->truth.index >= 0 );
    });
const Cut kIsCosmic = ( !kIsNuSlice );
const Cut kIsNuMuCC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.iscc && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ) );
    });
const Cut kIsNuOther = ( kIsNuSlice && !kIsNuMuCC );

const Cut kNoCut_numucc = kNoCut && kIsNuMuCC;
const Cut kNoCut_cosmic = kNoCut && kIsCosmic;
const Cut kNoCut_othernu = kNoCut && kIsNuOther;
const Cut kPreSelection_numucc = kSlcNuScoreCut && kInFV && kIsNuMuCC;
const Cut kPreSelection_cosmic = kSlcNuScoreCut && kInFV && kIsCosmic;
const Cut kPreSelection_othernu = kSlcNuScoreCut && kInFV && kIsNuOther;
const Cut kFlashMatching_numucc = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kIsNuMuCC;
const Cut kFlashMatching_cosmic = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kIsCosmic;
const Cut kFlashMatching_othernu = kSlcNuScoreCut && kInFV &&kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kIsNuOther;
const Cut kMuonTrk_numucc = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kIsNuMuCC;
const Cut kMuonTrk_cosmic = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kIsCosmic;
const Cut kMuonTrk_othernu = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kIsNuOther;
const Cut kALLCuts_numucc = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut && kIsNuMuCC;
const Cut kALLCuts_cosmic = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut && kIsCosmic;
const Cut kALLCuts_othernu = kSlcNuScoreCut && kInFV && kSlcFlashMatchTimeCut && kSlcFlashMatchScoreCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut && kIsNuOther;

// We are making the Spectra defined above
std::vector<SelDef> sels =
  {{"nocut_numucc", "All Slices",  kNoCut_numucc, kRed},
   {"nocut_cosmic", "All Slices",  kNoCut_cosmic, kRed+1},
   {"nocut_othernu", "All Slices",  kNoCut_othernu, kRed+2},
   {"PreSel_numucc",  "PreSel",  kPreSelection_numucc,  kOrange+7},
   {"PreSel_cosmic",  "PreSel",  kPreSelection_cosmic,  kOrange+6},
   {"PreSel_othernu",  "PreSel",  kPreSelection_othernu,  kOrange+5},
   {"FMatch_numucc",  "PreSel+FlashM",  kFlashMatching_numucc,  kMagenta-4},
   {"FMatch_cosmic",  "PreSel+FlashM",  kFlashMatching_cosmic,  kMagenta-3},
   {"FMatch_othernu",  "PreSel+FlashM",  kFlashMatching_othernu,  kMagenta-2},
   {"Mtrk_numucc",  "PreSel+FlashM+MuonID",  kMuonTrk_numucc,  kCyan-7},
   {"Mtrk_cosmic",  "PreSel+FlashM+MuonID",  kMuonTrk_cosmic,  kCyan-6},
   {"Mtrk_othernu",  "PreSel+FlashM+MuonID",  kMuonTrk_othernu,  kCyan-5},
   {"all_numucc",  "PreSel+FlashM+MuonID+CRT",  kALLCuts_numucc,  kGreen+2},
   {"all_cosmic",  "PreSel+FlashM+MuonID+CRT",  kALLCuts_cosmic,  kGreen+3},
   {"all_othernu",  "PreSel+FlashM+MuonID+CRT",  kALLCuts_othernu,  kGreen+4},
  };

// If you wanted to  add a new cut, define an
// expression for kYourCut (see examples in SBNAna/Cuts)
// and then add the following lines to the sels vector:
/* {"all_yourcut", "All Slices",  kYourCut,            kBlack}, */
/* {"sig_yourcut", "True NumuCC", (kSig && kYourCut),  kRed+1}, */
/* {"bkg_yourcut", "Not NumuCC",  (!kSig && kYourCut), kAzure+2}, */
