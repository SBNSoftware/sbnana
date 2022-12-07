#include <cmath>
#include <iostream>
#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Vars/NumuVarsSBND202106.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana{

  const Cut kInFV([](const caf::SRSliceProxy* slc)
                        {
                          return ( slc->vertex.x > (-199.15 + 10) && slc->vertex.x < (199.15 - 10) && slc->vertex.y > (-200. + 10) && slc->vertex.y < (200. - 10) && slc->vertex.z > (0.0 + 10) && slc->vertex.z < (500. - 50));
                        });

  const Cut kSlcFlashMatchTimeCut([](const caf::SRSliceProxy* slc)
                        {
                          if(std::isnan(slc->fmatch.time)) return false;

                          bool InBeam = (slc->fmatch.time > 0. && slc->fmatch.time < 1.800);
                          return ( InBeam );
                        });

  const Cut kSlcFlashMatchScoreCut([](const caf::SRSliceProxy* slc)
                        {
                          return ( slc->fmatch.score < 7 );
                        });

  //const Cut kHasPrimaryMuonTrk = kPrimaryMuonTrkIdx != -1;
  const Cut kHasPrimaryMuonTrk([](const caf::SRSliceProxy* slc)
                        {
                          int ptrkid = kPrimaryMuonTrkIdx(slc);
                          double ptrkrecop = kPrimaryMuonTrkP(slc);
                          //return ( kPrimaryMuonTrkIdx != -1 && !isnan(kPrimaryMuonTrkP) && kPrimaryMuonTrkP > 0. && kPrimaryMuonTrkP < 7.5);
                          return ( ptrkid != -1 && !isnan(ptrkrecop) && ptrkrecop > 0. && ptrkrecop < 7.5);
                        });


  //CRT Matching
  const Cut kCRTTrackAngleCut([](const caf::SRSliceProxy* slc)
                        {
                          int ptrkid = kPrimaryMuonTrkIdx(slc);
			  const caf::SRTrackProxy& ptrk = slc->reco.pfp[ptrkid].trk;
                          return ( isnan(ptrk.crttrack.angle) || ptrk.crttrack.angle > 0.05);
                        });

  const Cut kCRTHitDistanceCut([](const caf::SRSliceProxy* slc)
                        {
                          int ptrkid = kPrimaryMuonTrkIdx(slc);
			  const caf::SRCRTHitMatchProxy& crthit = slc->reco.pfp[ptrkid].trk.crthit;

                          // No candidate CRT hit matches
                          if(std::isnan(crthit.hit.time) || std::isnan(crthit.distance)) return true;

			  //Beam time window used in SBNSoftware/sbncode/blob/develop/sbncode/NuMuSelection/jupyter-ana/selection.py
                          return crthit.distance > 5 || (crthit.hit.time > 0. && crthit.hit.time < 1.800);
                        });

}
