#include "sbnana/SBNAna/Vars/NumuVarsSBND202106.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"
#include <cassert>

namespace ana
{


//  const Var kIsTrackAtSlc([](const caf::SRSliceProxy *slc, const FidVol& vol)
//			       {
//  				 double dist = sqrt((slc->reco.trk.start.x - slc->vertex.x)^2 + (slc->reco.trk.start.y - slc->vertex.y)^2 + (slc->reco.trk.start.z - slc->vertex.z)^2);
//  				 return (dist < 10); 
//  				
//			       });

//  const Var kIsTrackContained([](const caf::SRSliceProxy *slc, const FidVol& vol)
//			       {
//  				 return (vol.xmin < slc->reco.trk.end.x && slc->reco.trk.end.x < vol.xmax &&
//          				 vol.ymin < slc->reco.trk.end.y && slc->reco.trk.end.y < vol.ymax &&
//         				 vol.zmin < slc->reco.trk.end.z && slc->reco.trk.end.z < vol.zmax);
//			       });

  const Var kPrimaryMuonTrkIdx([](const caf::SRSliceProxy *slc, const FidVol& vol)
			       {
   				 double longest = -1;
   				 int best_idx = -1;

   				 for(unsigned int trkidx = 0; trkidx < slc->reco.trk.size(); ++trkidx){
     				   const caf::SRTrackProxy& trk = slc->reco.trk[trkidx];
				   
				   //atslc
  				   double dist = sqrt((pow(trk.start.x - slc->vertex.x,2)) + (pow(trk.start.y - slc->vertex.y,2)) + (pow(trk.start.z  - slc->vertex.z,2)));
  				   bool atslc = dist < 10;
     				   if(!atslc || !trk.parent_is_primary) continue;

  				   //contained
  				   bool contained = ( vol.xmin < trk.end.x && trk.end.x < vol.xmax && vol.ymin < trk.end.y && trk.end.y < vol.ymax && vol.zmin < trk.end.z && trk.end.z < vol.zmax);

				   //bestplane



     				   bool maybe_muon_exiting = !contained && trk.len > 100;  
				   bool maybe_muon_contained = contained && trk.bestplane.chi2_proton > 60 && trk.bestplane.chi2_muon < 30 && trk.len > 50;

     				   if(!maybe_muon_contained && !maybe_muon_exiting) continue;

     				   if(trk.len > longest){
       				     longest = trk.len;
			             best_idx = trkidx;
     				   }	
   				 }
   			         return best_idx;
			       });
}
