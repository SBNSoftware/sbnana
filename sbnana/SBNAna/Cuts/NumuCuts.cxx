#include <iostream>
#include "sbnana/SBNAna/Cuts/NumuCuts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana{

  const Cut kNumuBasicQual([](const caf::SRSliceProxy* slc)
			   { // 
			     bool hastrk = ( slc->reco.npfp > 0 );
			     if (!hastrk) return hastrk;

			     int muIdx = kPrimMuonIdx(slc);
			     if  (muIdx < 0) { return false; }
 			     double len = slc->reco.pfp[muIdx].trk.len;
			     return (len > 0 );
			   });

  const Cut kHasFlashMatch([](const caf::SRSliceProxy* slc)
			{ 
			  return ( slc->fmatch.present );
			});

  const Cut kFlashMatchScore([](const caf::SRSliceProxy* slc)
			{ 
			  return ( slc->fmatch.time > 0 && 
				   slc->fmatch.score > 6   );
			});

  const Cut kFlashMatchNumuCut = kHasFlashMatch  && kFlashMatchScore;

  
  const Cut kNumuTrkLen([](const caf::SRSliceProxy* slc)
			{ // TO DO: Prim pandora trk len primary_track.length > fConfig.TrackLength
			  bool hastrk = ( slc->reco.npfp > 0 );
			  if (!hastrk) return hastrk;
			  
			  int muIdx = kPrimMuonIdx(slc);
			  if (muIdx < 0) { return false; }
			  double len = slc->reco.pfp[muIdx].trk.len;
			  return (len > 50 );
			});
  
  // const Cut kNumuContND([](const caf::SRSliceProxy* slc)
  // 			{ // TO DO: Trk end-point containment
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			  if (!hastrk) return hastrk;
  // 			  // CosmicContain from sbnanalysis OscRecoPostprocess
  // 			  //  ytop: 40	 ybottom: 15  zfront: 15   zback: 15
  // 			  int muIdx = kPrimMuonIdx(slc);
  //                      if (muIdx < 0) { return false; }
  // 			  return ( slc->reco.pfp[muIdx].trk.start.x > ( kNDTopX -15 ) &&
  // 				   slc->reco.pfp[muIdx].trk.start.x < 15 &&
  // 				   slc->reco.pfp[muIdx].trk.start.y > ( kNDTopY -40 ) &&
  // 				   slc->reco.pfp[muIdx].trk.start.y < 15 &&
  // 				   slc->reco.pfp[muIdx].trk.start.z > ( kNDBack -15 ) &&
  // 				   slc->reco.pfp[muIdx].trk.start.z < 15 &&
  // 				   slc->reco.pfp[muIdx].trk.end.x > ( kNDTopX -15 ) &&
  // 				   slc->reco.pfp[muIdx].trk.end.x < 15 &&
  // 				   slc->reco.pfp[muIdx].trk.end.y > ( kNDTopY -40 ) &&
  // 				   slc->reco.pfp[muIdx].trk.end.y < 15 &&
  // 				   slc->reco.pfp[muIdx].trk.end.z > ( kNDBack -15 ) &&
  // 				   slc->reco.pfp[muIdx].trk.end.z < 15 && )
  // 			});
  
  // const Cut kNumuContFD([](const caf::SRSliceProxy* slc)
  // 			{ // TO DO: Trk end-point containment
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  // const Cut kNumuVtxFidND([](const caf::SRSliceProxy* slc)
  // 			  { // TO DO: Vtx Containment
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			    return ( hastrk );

  // 			  });
  
  // const Cut kNumuVtxFidFD([](const caf::SRSliceProxy* slc)
  // 			  { // TO DO: Vtx Containment
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			    return ( hastrk );
  // 			  });  
    
  // const Cut kNumuMCSLen([](const caf::SRSliceProxy* slc)
  // 			{ // TO DO: Prim trk MCS length -  primary_track.mcs_momentum < 7. /*garbage value*/&& 			  (fConfig.MCSTrackLength <0. || 			   primary_track.length > fConfig.MCSTrackLength
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  
  // const Cut kNumuCRTHit([](const caf::SRSliceProxy* slc)
  // 			{ // TO DO
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  // const Cut kNumuCRTTrk([](const caf::SRSliceProxy* slc)
  // 			{ // TO DO
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});
  
  // const Cut kNumuCRTAct([](const caf::SRSliceProxy* slc)
  // 			{ // TO DO
  // 			  bool hastrk = ( slc->reco.ntrk > 0 );
  // 			  return ( hastrk );
  // 			});  

}
