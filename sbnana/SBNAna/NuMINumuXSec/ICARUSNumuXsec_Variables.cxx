#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  // SpillVar
  const SpillVar spillvarCountSpill([](const caf::SRSpillProxy *sr) ->int {
    return 0.;
  });
  // - Test
  const SpillVar spillvarTest([](const caf::SRSpillProxy *sr) ->int {
    double ret = 0.;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      //std::cout << "[spillvarTest] Slice index = " << i << ", slc.reco.pfp.size() = " << slc.reco.pfp.size() << std::endl;
      int nTrk=0, nShw=0;
      for(std::size_t ip(0); ip < slc.reco.pfp.size(); ++ip){
        bool IsTrack = slc.reco.pfp.at(ip).trackScore > 0.5;
        if(IsTrack) nTrk++;
        else nShw++;
      }
      //std::cout << "[spillvarTest] npfp = " << slc.reco.pfp.size() << ", (Trk, Shw) = (" << nTrk << ", " << nShw << ")" << std::endl;
      std::cout << nTrk << "\t" << nShw << std::endl;
      ret = (nTrk+nShw);
    }
    return ret;
  });
  // - PMT-CRT matching
  const SpillMultiVar spillvarOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      rets.push_back(opflash.firsttime);
    }
    return rets;
  });
  const SpillMultiVar spillvarValidOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      if(opflash.onbeamtime) rets.push_back(opflash.firsttime); // TODO I'm using a hacked version of onbeamtime..
    }
    return rets;
  });
  const SpillMultiVar spillvarInTimeOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> validTimes = spillvarValidOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : validTimes){
      if( cpmt.IsInTime(opt) ) rets.push_back(opt);
    }
    return rets;
  });
  const SpillMultiVar spillvarCRTPMTTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : intimeTimes){
      int crtHitIdx = cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
      if(crtHitIdx>=0){
        rets.push_back( sr->crt_hits.at(crtHitIdx).t1 - opt );
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarCRTPMTMatchingID([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    if(intimeTimes.size()>0){
      for(const auto& opt : intimeTimes){
        rets.push_back( cpmt.GetMatchID(opt, sr->crt_hits) );
      }
    }
    return rets;
  });

  const SpillVar spillvarCRTPMTMatchingEventID([](const caf::SRSpillProxy *sr)
  {
/*
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    if(intimeTimes.size()>0){
      int CountNoMatching = 0;
      int CountEnteringTrack = 0;
      int CountExitingTrack = 0;
      int CountUnknown = 0
      for(const auto& opt : intimeTimes){
        int this_ID = cpmt.GetMatchID(opt, sr->crt_hits);
        if(this_ID==0) CountNoMatching++;
        if(this_ID==1 || this_ID==2 || this_ID==3 || this_ID==6 || this_ID==7) CountEnteringTrack++;
        if(this_ID==4 || this_ID==5) CountExitingTrack++;
        if(this_ID==8) CountUnknown++;
      }
    }
    else{
      return 9;
    }
*/
    return 0;
  });

  // Var
  // - Flash matching
  const Var varFMScore([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.score)) return -2.;
    else if(slc->fmatch.score<0) return -1.;
    else return slc->fmatch.score;
  });
  const Var varFMTime([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.time)) return -62.;
    else if(slc->fmatch.time<-50) return -61.;
    else return slc->fmatch.time;
  });
  // - Longest track
  //   - index
  const Var varLongestTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int ret(-1);
    double lmax(-999.);

    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){

      bool IsTrack = slc->reco.pfp.at(i).trackScore > 0.5;
      if(!IsTrack) continue;

      const auto& trk = slc->reco.pfp.at(i).trk;

      bool pass = false;
      if( fv_track.isContained(trk.end.x, trk.end.y, trk.end.z) ){
        pass = trk.len>50.;
      }
      else{
        pass = trk.len>100.;
      }
      if(!pass) continue;

      if(trk.len>lmax){
        lmax = trk.len;
        ret = i;
      }
    }

    return ret;
  });
  //   - length
  const Var varLongestTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.len;
    }
    else{
      return -999.;
    }
  });
  //   - direction
  const Var varLongestTrackDirectionX([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.dir.x;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackDirectionY([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.dir.y;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackDirectionZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.dir.z;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackDirectionXZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return sqrt(trk.dir.x*trk.dir.x+trk.dir.z*trk.dir.z);
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackForceDownDirectionX([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.x*flip;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackForceDownDirectionY([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.y*flip;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackForceDownDirectionZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.z*flip;
    }
    else{
      return -999.;
    }
  });
}
