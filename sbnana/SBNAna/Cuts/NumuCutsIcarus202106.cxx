#include <cmath>
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  const Cut kIsNuSlice = ( kTruthIndex >= 0.f );
  
  const Cut kIsCosmic = ( !kIsNuSlice );
  
  const Cut kIsNuMuCC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.iscc && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ) );
    });
  
  const Cut kIsNuOther = ( kIsNuSlice && !kIsNuMuCC );

  const Cut kCryo0([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->vertex.x) && slc->vertex.x < 0 );
    });

  const Cut kTFiducial([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->truth.position.x) &&
	       ( ( slc->truth.position.x < -71.1 - 25 && slc->truth.position.x > -369.33 + 25 ) ||
		 ( slc->truth.position.x > 71.1 + 25 && slc->truth.position.x < 369.33 - 25 ) ) &&
	       !isnan(slc->truth.position.y) &&
	       ( slc->truth.position.y > -181.7 + 25 && slc->truth.position.y < 134.8 - 25 ) &&
	       !isnan(slc->truth.position.z) &&
	       ( slc->truth.position.z > -895.95 + 30 && slc->truth.position.z < 895.95 - 50 ) );
    });

  const Cut kRFiducial([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->vertex.x) &&
	       ( ( slc->vertex.x < -71.1 - 25 && slc->vertex.x > -369.33 + 25 ) ||
		 ( slc->vertex.x > 71.1 + 25 && slc->vertex.x < 369.33 - 25 ) ) &&
	       !isnan(slc->vertex.y) &&
	       ( slc->vertex.y > -181.7 + 25 && slc->vertex.y < 134.8 - 25 ) &&
	       !isnan(slc->vertex.z) &&
	       ( slc->vertex.z > -895.95 + 30 && slc->vertex.z < 895.95 - 50 ) );
    });

  const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
      return !slc->is_clear_cosmic;
    });

  const Cut kNuScore([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->nu_score) && slc->nu_score > 0.4 );
    });

  const Cut kFMScore([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 7.0 );
    });

  const Cut kPTrack([](const caf::SRSliceProxy* slc) {
      return ( kPTrackInd(slc) >= 0 );
    });

  const Cut kPTrackContained([](const caf::SRSliceProxy* slc) {
      int Ind = kPTrackInd(slc);
      bool Contained(false);
      if ( Ind >= 0 ) 
      {
	auto const& trk = slc->reco.pfp.at(Ind).trk;
	Contained = ( !isnan(trk.end.x) &&
		      ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
		      !isnan(trk.end.y) &&
		      ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		      !isnan(trk.end.z) &&
		      ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
      }
      return Contained;
    });

  const Cut kPTrackExiting([](const caf::SRSliceProxy* slc) {
      int Ind = kPTrackInd(slc);
      bool Exiting(false);
      if ( Ind >= 0 ) 
      {
	auto const& trk = slc->reco.pfp.at(Ind).trk;
	Exiting = !( !isnan(trk.end.x) &&
		     ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
		     !isnan(trk.end.y) &&
		     ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		     !isnan(trk.end.z) &&
		     ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
      }
      return Exiting;
    });

  // "Successive" cuts for the selection. These are meant to successively build upon
  // each other, which fits the flow of the selection.
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
  
  // The full selection cut using only reco information.
  const Cut kNuMuCC_FullSelection = ( kCryo0 && kRFiducial && kNotClearCosmic && kNuScore && kFMScore && kPTrack );

  const SpillCut kLongTrack([](const caf::SRSpillProxy* sr){
      bool ProgCut = false;
      bool NuIsNuMuCC, IsMuon, Contained;
      for (auto const& nu: sr->mc.nu)
      {
	NuIsNuMuCC = nu.iscc && ( nu.pdg == 14 || nu.pdg == -14 );
	for (auto const& prim: nu.prim)
	{
	  IsMuon = ( prim.pdg == 13 || prim.pdg == -13 );
	  Contained = ( prim.contained == 1 );
	  ProgCut = ProgCut || ( NuIsNuMuCC && IsMuon && ( (Contained && prim.length > 50.0) || (!Contained && prim.length > 100.0) ) );
	}
      }
      return ProgCut;
    });
}
