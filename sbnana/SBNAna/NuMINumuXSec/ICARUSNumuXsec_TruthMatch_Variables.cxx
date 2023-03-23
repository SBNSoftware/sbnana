#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TruthMatch{

  // For a given true muon (truth_index), find a reco track whose best-matched is this particle
  const Var TruthMuonIndex([](const caf::SRSliceProxy* slc) -> int {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg) == 13 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          truth_idx = i;
        }
      }
    }
    return truth_idx;
  });
  const Var TruthMuonLength([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonIndex(slc);
    if(truth_idx>=0){
      if(isnan(slc->truth.prim.at(truth_idx).length)) return -999.;
      else return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonKE([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).genE - M_MUON;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int truth_idx = TruthMuonIndex(slc);
    return GetMatchedRecoTrackIndex(slc, truth_idx);
  });
  const Var TruthMuonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    // -1 : No track found
    // 0 : Exiting
    // 1 : Contained
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var TruthMuonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
    }
    else{
      return 999999;
    }
  });

  const Var TruthMuonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      if(isnan(Chi2Muon)) return -999.;
      else return Chi2Muon;
    }
    else{
      return 999999;
    }
  });
  const Var TruthMuonMatchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      if(isnan(Chi2Pion)) return -999.;
      else return Chi2Pion;
    }
    else{
      return 999999;
    }
  });

  // For a given true proton (truth_index), find a reco track whose best-matched is this particle
  const Var TruthProtonIndex([](const caf::SRSliceProxy* slc) -> int {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( slc->truth.prim.at(i).pdg == 2212 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          truth_idx = i;
        }
      }
    }
    return truth_idx;
  });
  const Var TruthProtonLength([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthProtonIndex(slc);
    if(truth_idx>=0){
      if(isnan(slc->truth.prim.at(truth_idx).length)) return -999.;
      else return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }
  });
  const Var TruthProtonKE([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthProtonIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).genE - M_PROTON;
    }
    else{
      return -999.;
    }
  });
  const Var TruthProtonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int truth_idx = TruthProtonIndex(slc);
    return GetMatchedRecoTrackIndex(slc, truth_idx);
  });
  const Var TruthProtonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    // -1 : No track found
    // 0 : Exiting
    // 1 : Contained
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var TruthProtonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
    }
    else{
      return 999999;
    }
  });

  const Var TruthProtonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      if(isnan(Chi2Muon)) return -999.;
      else return Chi2Muon;
    }
    else{
      return 999999;
    }
  });
  const Var TruthProtonMatchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      if(isnan(Chi2Pion)) return -999.;
      else return Chi2Pion;
    }
    else{
      return 999999;
    }
  });

  // For a given true charged pion (truth_index), find a reco track whose best-matched is this particle
  const Var TruthChargedPionIndex([](const caf::SRSliceProxy* slc) -> int {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg) == 211){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          truth_idx = i;
        }
      }
    }
    return truth_idx;
  });
  const Var TruthChargedPionLength([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      if(isnan(slc->truth.prim.at(truth_idx).length)) return -999.;
      else return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionKE([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).genE - M_CHARGEDPION;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx<0) return -999.;
    return GetMatchedRecoTrackIndex(slc, truth_idx);
  });
  const Var TruthChargedPionMatchedTrackEndProcess([](const caf::SRSliceProxy* slc) -> int {
    int trackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(trackIndex>=0){
      const auto& trk = slc->reco.pfp.at(trackIndex).trk;
      return trk.truth.p.end_process;
    }
    else{
      return -999.;
    }
  });


} // end namespace TruthMatch

} // end namespace ICARUSNumuXsec
