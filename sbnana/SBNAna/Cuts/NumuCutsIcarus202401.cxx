#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202401.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <iostream>

namespace ana {

static bool Icarus202401contained(const caf::SRTrackProxy& trk)
{
  const double x = trk.end.x, y = trk.end.y, z = trk.end.z;
  if(std::isnan(x) || std::isnan(y) || std::isnan(z)) return false;

  bool x_contained_E = x <  -61.94 - 5 && x > -358.49 + 5;
  bool x_contained_W = x >   61.94 + 5 && x <  358.49 - 5;
  bool y_contained   = y > -181.86 + 5 && y <  134.96 - 5;
  bool z_contained   = z > -894.95 + 5 && z <  894.95 - 5;

  // Check triangles at corners of the detector
  bool in_corner = y <  1.732007 * z - 1687.5114 
                || y > -1.732007 * z + 1640.6114
                || y >  1.732007 * z + 1640.6114 
                || y < -1.732007 * z - 1687.5114;

  return (x_contained_E || x_contained_W) && y_contained && z_contained 
         && !in_corner;
}

const Cut kIcarus202401BaryFMCut([](const caf::SRSliceProxy *slc) {
    return !std::isnan(slc->barycenterFM.deltaZ_Trigger) && 
           slc->barycenterFM.deltaZ_Trigger >= 0 && 
           slc->barycenterFM.deltaZ_Trigger < 100;
});

const SpillCut kIcarus202401CRTPMTVeto([](const caf::SRSpillProxy* spill){
    double min_time = 0, max_time = 1.6;
    if(!spill->hdr.ismc){min_time = -0.4; max_time = 1.5;}
    for(const auto& match: spill->crtpmt_matches) {
        if(match.flashGateTime > min_time && 
           match.flashGateTime < max_time && 
           match.flashClassification == 0) {
            return true;
        } 
    }
    return false;
});

const Cut kIcarus202401FoundMuon = kIcarus202401MuonIdx >= 0;
const Cut kIcarus202401RecoFiducial([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->vertex.x) &&
	       ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
		 ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
	       !isnan(slc->vertex.y) &&
	       ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
	       !isnan(slc->vertex.z) &&
	       ( slc->vertex.z > -894.95 + 30 && slc->vertex.z < 894.95 - 50 ) &&
      !(slc->vertex.x > 210 && slc->vertex.y > 60 && slc->vertex.z > 290 && slc->vertex.z < 390));
    });

const Cut kIcarus202401NumuSelection = kIcarus202401RecoFiducial && kIcarus202401BaryFMCut && kIcarus202401FoundMuon;

const Cut kIcarus202401NoPion = kIcarus202401NumPions == 0;
const Cut kIcarus202401ContainedHadrons([](const caf::SRSliceProxy* slc){
  auto idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    //if (pfp.trackScore < 0.5) { continue; }

    if(std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x) || std::isnan(pfp.trk.len)) continue;
    auto const& trk = pfp.trk;
    //if(pfp.id != muID && pfp.parent_is_primary)
    if(pfp.id != muID)
      if(!Icarus202401contained(trk)) return false;
  }
  return true;
});
const Cut kIcarus202401ContainedMuon([](const caf::SRSliceProxy* slc){
  return kIcarus202401FoundMuon(slc) && Icarus202401contained(slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk);
});
const Cut kIcarus202401ContainedMuonAndHadrons = kIcarus202401ContainedMuon && kIcarus202401ContainedHadrons;
const Cut kIcarus202401ContainedNumuCCInclusive = kIcarus202401NumuSelection && kIcarus202401ContainedHadrons && kIcarus202401ContainedMuon;
const Cut kIcarus202401NumuCC0pi = kIcarus202401ContainedNumuCCInclusive && kIcarus202401NoPion;
const Cut kIcarus202401Contained1muNp = kIcarus202401NumuCC0pi && kIcarus202401NumProtons > 0 && kIcarus202401NumShowers == 0;
const Cut kIcarus202401Contained1mu1p = kIcarus202401Contained1muNp && kIcarus202401NumProtons == 1;
}
