// SBNAna includes.
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

// Custom includes.
#include "TestVars.h"

using namespace ana;

////////////////////////////////////////////////////////////////////////////////////////////////////
// Define the cuts for the selection.
////////////////////////////////////////////////////////////////////////////////////////////////////

const Cut kIsNu = ( kTruthIndex >= 0.f );

const Cut kIsCosmic = ( !kIsNu );

const Cut kIsNuMuCC([](const caf::SRSliceProxy* slc) {
    return ( kIsNu(slc) && slc->truth.iscc && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ) );
  });

const Cut kIsNuOther = ( kIsNu && !kIsNuMuCC );

const Cut kCryo0([](const caf::SRSliceProxy* slc) {
    return ( !std::isnan(slc->vertex.x) && slc->vertex.x < 0 );
  });

const Cut kTFiducial([](const caf::SRSliceProxy* slc) {
    return ( !std::isnan(slc->truth.position.x) &&
	     ( ( slc->truth.position.x < -71.1 - 25 && slc->truth.position.x > -369.33 + 25 ) ||
	       ( slc->truth.position.x > 71.1 + 25 && slc->truth.position.x < 369.33 - 25 ) ) &&
	     !std::isnan(slc->truth.position.y) &&
	     ( slc->truth.position.y > -181.7 + 25 && slc->truth.position.y < 134.8 - 25 ) &&
	     !std::isnan(slc->truth.position.z) &&
	     ( slc->truth.position.z > -895.95 + 30 && slc->truth.position.z < 895.95 - 50 ) );
  });

const Cut kRFiducial([](const caf::SRSliceProxy* slc) {
    return ( !std::isnan(slc->vertex.x) &&
	     ( ( slc->vertex.x < -71.1 - 25 && slc->vertex.x > -369.33 + 25 ) ||
	       ( slc->vertex.x > 71.1 + 25 && slc->vertex.x < 369.33 - 25 ) ) &&
	     !std::isnan(slc->vertex.y) &&
	     ( slc->vertex.y > -181.7 + 25 && slc->vertex.y < 134.8 - 25 ) &&
	     !std::isnan(slc->vertex.z) &&
	     ( slc->vertex.z > -895.95 + 30 && slc->vertex.z < 895.95 - 50 ) );
  });

const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !slc->is_clear_cosmic;
  });

const Cut kNuScore([](const caf::SRSliceProxy* slc) {
    return ( !std::isnan(slc->nu_score) && slc->nu_score > 0.4 );
  });

const Cut kFMScore([](const caf::SRSliceProxy* slc) {
    return ( !std::isnan(slc->fmatch.score) && slc->fmatch.score < 7.0 );
  });

const Cut kPTrack([](const caf::SRSliceProxy* slc) {
    // The (dis)qualification of a slice is based upon the track level information.
    float Atslc, Chi2Proton, Chi2Muon;
    std::vector<unsigned int> Candidates;
    bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained, ProgCut(false);

    for (auto const& trk : slc->reco.trk)
    {
      // First we calculate the distance of each track to the slice vertex.
      Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
		    pow( slc->vertex.y - trk.start.y, 2.0 ) +
		    pow( slc->vertex.z - trk.start.z, 2.0 ) );
      // We require that the distance of the track from the slice is less than
      // 10 cm and that the parent of the track has been marked as the primary.
      AtSlice = ( Atslc < 10.0 && trk.parent_is_primary);
      
      if (trk.bestplane == 0)
      {
	Chi2Proton = trk.chi2pid0.chi2_proton;
	Chi2Muon = trk.chi2pid0.chi2_muon;
      }
      else if (trk.bestplane == 1)
      {
	Chi2Proton = trk.chi2pid1.chi2_proton;
	Chi2Muon = trk.chi2pid1.chi2_muon;
      }
      else
      {
	Chi2Proton = trk.chi2pid2.chi2_proton;
	Chi2Muon = trk.chi2pid2.chi2_muon;
      }

      Contained = ( !std::isnan(trk.end.x) &&
		    ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&//
		      // || ( trk.end.x > 71.1 + 25 && trk.end.x < 369.33 - 25 ) ) &&
		    !std::isnan(trk.end.y) &&
		    ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		    !std::isnan(trk.end.z) &&
		    ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
      MaybeMuonExiting = ( !Contained && trk.len > 100);
      MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
      ProgCut = ProgCut || ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) );
    }
    return ProgCut;
  });

const Cut kNuMuCC_Cryo0 = ( kIsNuMuCC && kCryo0 );
const Cut kCosmic_Cryo0 = ( kIsCosmic && kCryo0 );
const Cut kNuOther_Cryo0 = ( kIsNuOther && kCryo0 );

const Cut kNuMuCC_TFiducial = ( kNuMuCC_Cryo0 && kTFiducial );
const Cut kNuMuCC_RFiducial = ( kNuMuCC_Cryo0 && kRFiducial );
const Cut kCosmic_TFiducial = ( kCosmic_Cryo0 && kTFiducial );
const Cut kCosmic_RFiducial = ( kCosmic_Cryo0 && kRFiducial );
const Cut kNuOther_TFiducial = ( kNuOther_Cryo0 && kTFiducial );
const Cut kNuOther_RFiducial = ( kNuOther_Cryo0 && kRFiducial );

const Cut kNuMuCC_ClearCos = ( kNuMuCC_TFiducial && kNotClearCosmic );
const Cut kCosmic_ClearCos = ( kCosmic_RFiducial && kNotClearCosmic );
const Cut kNuOther_ClearCos = ( kNuOther_TFiducial && kNotClearCosmic );

const Cut kNuMuCC_NuScore = ( kNuMuCC_ClearCos && kNuScore );
const Cut kCosmic_NuScore = ( kCosmic_ClearCos && kNuScore );
const Cut kNuOther_NuScore = ( kNuOther_ClearCos && kNuScore );

const Cut kNuMuCC_FMScore = ( kNuMuCC_NuScore && kFMScore );
const Cut kCosmic_FMScore = ( kCosmic_NuScore && kFMScore );
const Cut kNuOther_FMScore = ( kNuOther_NuScore && kFMScore );

const Cut kNuMuCC_PTrack = ( kNuMuCC_FMScore && kPTrack );
const Cut kCosmic_PTrack = ( kCosmic_FMScore && kPTrack );
const Cut kNuOther_PTrack = ( kNuOther_FMScore && kPTrack );
