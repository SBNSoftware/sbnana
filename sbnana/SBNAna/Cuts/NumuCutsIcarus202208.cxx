#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202208.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202208.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana {

static bool Icarus202208contained(const caf::SRTrackProxy& trk)
{
  return ((trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5)
             || (trk.end.x > 61.94 + 5 && trk.end.x < 358.49 - 5))
           && trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5
           && trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5;
}

const Cut kIcarus202208FMTimeCut = kFMTimeVar > 0 && kFMTimeVar < 1.8;
const Cut kIcarus202208FMScoreCut = kFMScoreVar < 9;
const Cut kIcarus202208LongTrackDirCut = kCRLongestTrackDirY > -0.91;
const Cut kIcarus202208FoundMuon = kIcarus202208MuonIdx >= 0;
const Cut kIcarus202208RecoFiducial([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->vertex.x) &&
	       ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
		 ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
	       !isnan(slc->vertex.y) &&
	       ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
	       !isnan(slc->vertex.z) &&
	       ( slc->vertex.z > -894.95 + 30 && slc->vertex.z < 894.95 - 50 ) );
    });

const Cut kIcarus202208NumuSelection = kIcarus202208RecoFiducial && kIcarus202208FMScoreCut && kIcarus202208FMTimeCut && kIcarus202208LongTrackDirCut && kIcarus202208FoundMuon;

const Cut kIcarus202208NoPion = kIcarus202208NumPions == 0;
const Cut kIcarus202208ContainedHadrons([](const caf::SRSliceProxy* slc){
  auto idx = kIcarus202208MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    if (pfp.trackScore < 0.5) { continue; }
    auto const& trk = pfp.trk;
    if(pfp.id != muID && pfp.parent_is_primary)
      if(!Icarus202208contained(trk)) return false;
  }
  return true;
});
const Cut kIcarus202208ContainedMuon([](const caf::SRSliceProxy* slc){
  return kIcarus202208FoundMuon(slc) && Icarus202208contained(slc->reco.pfp.at(kIcarus202208MuonIdx(slc)).trk);
});
const Cut kIcarus202208ContainedMuonAndHadrons = kIcarus202208ContainedMuon && kIcarus202208ContainedHadrons;
const Cut kIcarus202208QELike = kIcarus202208NumuSelection && kIcarus202208NoPion && kIcarus202208ContainedHadrons;
const Cut kIcarus202208QELikeContainedMuon = kIcarus202208QELike && kIcarus202208ContainedMuon;
}
