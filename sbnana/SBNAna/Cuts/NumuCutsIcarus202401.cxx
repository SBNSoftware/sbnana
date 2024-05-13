#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202401.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana {

static bool Icarus202401contained(const caf::SRTrackProxy& trk)
{
  return ((trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5)
             || (trk.end.x > 61.94 + 5 && trk.end.x < 358.49 - 5))
           && trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5
           && trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5;
}

const SpillCut kIcarus202401CRTPMTVeto([](const caf::SRSpillProxy* spill){
    for(const auto& match: spill->crtpmt_matches) {
        if(match.flashGateTime > 0 && match.flashGateTime < 1.6 && match.flashClassification == 0)
            return true;
    }
    return false;
});

const SpillCut kIcarus202401CRTPMTVetoData([](const caf::SRSpillProxy* spill){
    for(const auto& match: spill->crtpmt_matches) {
        if(match.flashGateTime > -0.4 && match.flashGateTime < 1.6 && match.flashClassification == 0)
            return true;
    }
    return false;
});

const Cut kIcarus202401BaryFMCut([](const caf::SRSliceProxy *slc) {
    return !std::isnan(slc->barycenterFM.deltaZ_Trigger) && 
           slc->barycenterFM.deltaZ_Trigger >= 0 && 
           slc->barycenterFM.deltaZ_Trigger < 100;
});

const Cut kIcarus202401FoundMuon = kIcarus202401MuonIdx >= 0;
const Cut kIcarus202401RecoFiducial([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->vertex.x) &&
	       ( ( slc->vertex.x < -61.94 - 25 && slc->vertex.x > -358.49 + 25 ) ||
		 ( slc->vertex.x > 61.94 + 25 && slc->vertex.x < 358.49 - 25 ) ) &&
	       !isnan(slc->vertex.y) &&
	       ( slc->vertex.y > -181.86 + 25 && slc->vertex.y < 134.96 - 25 ) &&
	       !isnan(slc->vertex.z) &&
	       ( slc->vertex.z > -894.95 + 30 && slc->vertex.z < 894.95 - 50 ) );
    });

const Cut kIcarus202401NumuSelection = kIcarus202401RecoFiducial && kIcarus202401BaryFMCut && kIcarus202401FoundMuon;

const Cut kIcarus202401NoPion = kIcarus202401NumPions == 0;
const Cut kIcarus202401ContainedHadrons([](const caf::SRSliceProxy* slc){
  auto idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    if (pfp.trackScore < 0.5) { continue; }
    auto const& trk = pfp.trk;
    if(pfp.id != muID && pfp.parent_is_primary)
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
