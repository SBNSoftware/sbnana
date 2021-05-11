#include <iostream>
#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Vars/NumuVarsSBND202106.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana{

  const Cut kHasPrimaryMuonTrk = kPrimaryMuonTrkIdx != -1;

  //CRT Matching
  const Cut kCRTTrackAngleCut([](const caf::SRSliceProxy* slc)
                        {
                          int ptrkid = kPrimaryMuonTrkIdx(slc);
			  const caf::SRTrackProxy& ptrk = slc->reco.trk[ptrkid];
                          return ( ptrk.crttrack.angle > 0.05);
                        });

  const Cut kCRTHitDistanceCut([](const caf::SRSliceProxy* slc)
                        {
                          int ptrkid = kPrimaryMuonTrkIdx(slc);
			  const caf::SRTrackProxy& ptrk = slc->reco.trk[ptrkid];
			  //Beam time window used in SBNSoftware/sbncode/blob/develop/sbncode/NuMuSelection/jupyter-ana/selection.py
                          return ( ptrk.crthit.distance > 5 && ptrk.crthit.hit.time > 0. && ptrk.crthit.hit.time < 1.800 );
                        });

}
