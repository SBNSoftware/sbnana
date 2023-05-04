#include "sbnana/SBNAna/Vars/NumuVarsIcarus202208.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana {

// More general vars, i.e., not just numu vars.
// Probably should be moved to another file at some point.
const Var kCRLongestTrackDirY = SIMPLEVAR(nuid.crlongtrkdiry);
const Var kFMTimeVar          = SIMPLEVAR(fmatch.time);
const Var kFMScoreVar         = SIMPLEVAR(fmatch.score);

// Mostly just copied from NumuVarsIcarus202106.cxx
const Var kIcarus202208MuonIdx([](const caf::SRSliceProxy* slc) -> int {
      // The (dis)qualification of a slice is based upon the track level information.
      float Longest(0);
      int PTrackInd(-1);
      for (std::size_t i(0); i < slc->reco.npfp; ++i) {
        auto const& pfp = slc->reco.pfp.at(i);
        if (pfp.trackScore < 0.5) { continue; }
        auto const& trk = pfp.trk;
        if(trk.bestplane == -1) continue;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

        const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
        const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
  
        const bool Contained = ( !isnan(trk.end.x) &&
        (( trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5 ) ||
        ( trk.end.x > 61.94 + 5 && trk.end.x < +358.49 - 5 )) &&
        !isnan(trk.end.y) &&
        ( trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5 ) &&
        !isnan(trk.end.z) &&
        ( trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5 ) );
        const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest ) {
          Longest = trk.len;
          PTrackInd = i;
        }
      }
      return PTrackInd;
});

static bool Icarus202208_proton_cut(const caf::SRTrackProxy& trk)
{
  return trk.chi2pid[trk.bestplane].chi2_proton < 100;
}

const Var kIcarus202208NumPions([](const caf::SRSliceProxy* slc)
{
  int count = 0;
  auto idx = kIcarus202208MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    if (pfp.trackScore < 0.5) { continue; }
    auto const& trk = pfp.trk;
    if(pfp.id != muID && !Icarus202208_proton_cut(trk) && pfp.parent_is_primary)
      ++count;
  }
  return count;
});
}
