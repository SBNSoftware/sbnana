#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_ThreeTrack_Cuts.h"

using namespace ana;
using namespace std;
using namespace ICARUSNumuXsec::TwoTrack;

namespace ICARUSNumuXsec{

namespace ThreeTrack{

  const Cut Pass3T1PCut([](const caf::SRSliceProxy* slc) {

    // count stuff
    unsigned n_proton = 0;
    unsigned n_trk_proton = 0;

    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
      const auto& pfp = slc->reco.pfp.at(i_pfp);

      bool is_track = pfp.trackScore > 0.5;
      const auto &pid = pfp.trk.chi2pid[2];
      bool is_proton = (pid.pid_ndof > 0) && (pid.chi2_proton < 60) && (pid.chi2_muon > 35);
      n_trk_proton += (is_track | is_proton);
      n_proton += is_proton;
    }
    bool pass = n_proton >= 1 && n_trk_proton >= 3 && !(slc->is_clear_cosmic);
    return pass;

  });

  const Cut HasThreePrimaryTracks([](const caf::SRSliceProxy* slc) {
    return PrimaryTrackIndices(slc).size()>=3;
  });
  const Cut HasMuonTrack([](const caf::SRSliceProxy* slc) {
    return MuonTrackIndex(slc)>=0.;
  });

  const Cut HasProtonTrack([](const caf::SRSliceProxy* slc) {
    return ProtonTrackIndex(slc)>=0;
  });

  const Cut MuonContained([](const caf::SRSliceProxy* slc) {
    int longerMuonTrackIndex = MuonTrackIndex(slc);
    if(longerMuonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(longerMuonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      return isContained;
    }
    else{
      return false;
    }
  });
  const Cut MuonExiting([](const caf::SRSliceProxy* slc) {
    int longerMuonTrackIndex = MuonTrackIndex(slc);
    if(longerMuonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(longerMuonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      return !isContained;
    }
    else{
      return false;
    }
  });

} // end namespace ThreeTrack

} // end namespace ICARUSNumuXsec
