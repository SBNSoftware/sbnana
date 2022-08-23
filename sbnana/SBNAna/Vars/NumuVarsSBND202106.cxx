#include "sbnana/SBNAna/Vars/NumuVarsSBND202106.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
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

  const Var kPrimaryMuonTrkIdx([](const caf::SRSliceProxy *slc)
			       {
   				 double longest = -1;
   				 int best_idx = -1;
				 double dist = -1;
                                 bool atslc, contained, maybe_muon_exiting, maybe_muon_contained;
				 float chi2_proton, chi2_muon;

   				 for(unsigned int trkidx = 0; trkidx < slc->reco.pfp.size(); ++trkidx){
     				   const caf::SRTrackProxy& trk = slc->reco.pfp[trkidx].trk;
                                   if(trk.bestplane == -1 || slc->reco.pfp[trkidx].trackScore < 0.5) continue;
				   
				   //atslc
  				   dist = sqrt((pow(trk.start.x - slc->vertex.x,2)) + (pow(trk.start.y - slc->vertex.y,2)) + (pow(trk.start.z  - slc->vertex.z,2)));
  				   atslc = dist < 10;
     				   if(!atslc || !slc->reco.pfp[trkidx].parent_is_primary) continue;

  				   //contained
  				   //contained = ( vol.xmin < trk.end.x && trk.end.x < vol.xmax && vol.ymin < trk.end.y && trk.end.y < vol.ymax && vol.zmin < trk.end.z && trk.end.z < vol.zmax);
  				   //The following definition reproduce the numbers from SBNSoftware/sbncode/blob/develop/sbncode/NuMuSelection/jupyter-ana/selection.py
  				   //instead of using the defined Fiducial Volume check implemented in sbnana/SBNAna/Cuts/VolumeDefinitions.h 
  				   contained = ( (-199.15 + 10) < trk.end.x && trk.end.x < (199.15 - 10) && (-200. + 10) < trk.end.y && trk.end.y < (200. - 10) && (0.0 + 10) < trk.end.z && trk.end.z < (500. - 50));

				   //bestplane.chi2
                                   chi2_proton = trk.chi2pid[trk.bestplane].chi2_proton;
                                   chi2_muon = trk.chi2pid[trk.bestplane].chi2_muon;

     				   maybe_muon_exiting = !contained && trk.len > 100;  
				   maybe_muon_contained = contained && chi2_proton > 60 && chi2_muon < 30 && trk.len > 50;

     				   if(!maybe_muon_contained && !maybe_muon_exiting) continue;

     				   if(trk.len > longest){
       				     longest = trk.len;
			             best_idx = trkidx;
     				   }	
   				 }
   			         return best_idx;
			       });

  const Var kPrimaryMuonTrkP([](const caf::SRSliceProxy *slc)
			       {
  				 //bool contained; 
				 //int ptrkid = kPrimaryMuonTrkIdx(slc);

				 //const caf::SRTrackProxy& ptrk = slc->reco.trk[ptrkid];
				 //contained = ( (-199.15 + 10) < ptrk.end.x && ptrk.end.x < (199.15 - 10) && (-200. + 10) < ptrk.end.y && ptrk.end.y < (200. - 10) && (0.0 + 10) < ptrk.end.z && ptrk.end.z < (500. - 50));
				 //if(contained) return ptrk.rangeP.p_muon;
				 //else return ptrk.mcsP.fwdP_muon;                               

			         float recop(-5.f);
      				 bool contained(false);
      
			         if ( kPrimaryMuonTrkIdx(slc) >= 0 ){
			       	   auto const& ptrk = slc->reco.pfp.at(kPrimaryMuonTrkIdx(slc)).trk;
				   contained = ( (-199.15 + 10) < ptrk.end.x && ptrk.end.x < (199.15 - 10) && (-200. + 10) < ptrk.end.y && ptrk.end.y < (200. - 10) && (0.0 + 10) < ptrk.end.z && ptrk.end.z < (500. - 50));

			 	   if(contained) recop = ptrk.rangeP.p_muon;
				   else recop = ptrk.mcsP.fwdP_muon;
      				 }
      				 return recop;
			       });

  const Var kCRTTrkTime([](const caf::SRSliceProxy *slc) -> double
  			  {
			    float crttrktime(-5.f);
                            if ( kPrimaryMuonTrkIdx(slc) >= 0 ){
                              int ptrkid = kPrimaryMuonTrkIdx(slc);
                              const caf::SRTrackProxy& ptrk = slc->reco.pfp[ptrkid].trk;
			      crttrktime = ptrk.crthit.hit.time;
			    }
         		    return crttrktime;
       			  });

  const Var kCRTTrkAngle([](const caf::SRSliceProxy *slc) -> double
  			  {
			    float crttrkangle(-5.f);
                            if ( kPrimaryMuonTrkIdx(slc) >= 0 ){
                              int ptrkid = kPrimaryMuonTrkIdx(slc);
                              const caf::SRTrackProxy& ptrk = slc->reco.pfp[ptrkid].trk;
			      crttrkangle = ptrk.crttrack.angle;
			    }
         		    return crttrkangle;
       			  });

  const Var kCRTHitDist([](const caf::SRSliceProxy *slc) -> double
  			  {
			    float crttrkdist(-5.f);
                            if ( kPrimaryMuonTrkIdx(slc) >= 0 ){
                              int ptrkid = kPrimaryMuonTrkIdx(slc);
                              const caf::SRTrackProxy& ptrk = slc->reco.pfp[ptrkid].trk;
			      crttrkdist = ptrk.crthit.distance;
			    }
         		    return crttrkdist;
       			  });

  const Var kNuScore([](const caf::SRSliceProxy *slc) -> double
  			  {
         		    return slc->nu_score;
       			  });

  const Var kPrimaryMuonTrkLen([](const caf::SRSliceProxy *slc) -> double
  			  {
			    float ptrklen(-5.f);
                            if ( kPrimaryMuonTrkIdx(slc) >= 0 ){
                              int ptrkid = kPrimaryMuonTrkIdx(slc);
                              const caf::SRTrackProxy& ptrk = slc->reco.pfp[ptrkid].trk;
			      ptrklen = ptrk.len;
			    }
         		    return ptrklen;
       			  });

}
