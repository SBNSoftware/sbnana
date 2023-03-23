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
      int nTrk=0, nShw=0;
      const auto& slc = sr->slc.at(i);
      std::cout << "[spillvarTest] Slice index = " << i << ", slc.reco.pfp.size() = " << slc.reco.pfp.size() << std::endl;
      for(std::size_t ip(0); ip < slc.reco.pfp.size(); ++ip){
        bool IsTrack = slc.reco.pfp.at(ip).trackScore > 0.5;
        if(IsTrack) nTrk++;
        else nShw++;
      }
      ret = (nTrk+nShw);
      std::cout << "[spillvarTest] npfp = " << slc.reco.pfp.size() << ", (Trk, Shw) = (" << nTrk << ", " << nShw << ")" << std::endl;
    }
    return ret;
  });
  const SpillVar spillvarNTrack([](const caf::SRSpillProxy *sr) ->int {
    int nTrk=0;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      //std::cout << "[spillvarTest] Slice index = " << i << ", slc.reco.pfp.size() = " << slc.reco.pfp.size() << std::endl;
      for(std::size_t ip(0); ip < slc.reco.pfp.size(); ++ip){
        bool IsTrack = slc.reco.pfp.at(ip).trackScore > 0.5;
        if(IsTrack) nTrk++;
      }
    }
    return nTrk;
  });
  const SpillVar spillvarNShower([](const caf::SRSpillProxy *sr) ->int {
    int nShw=0;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      //std::cout << "[spillvarTest] Slice index = " << i << ", slc.reco.pfp.size() = " << slc.reco.pfp.size() << std::endl;
      for(std::size_t ip(0); ip < slc.reco.pfp.size(); ++ip){
        bool IsShower = slc.reco.pfp.at(ip).trackScore < 0.5;
        if(IsShower) nShw++;
      }
    }
    return nShw;
  });
  // - CRT Hit
  const SpillMultiVar spillvarSideCRTHitPe([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=40&&hit.plane<=49){
        rets.push_back(hit.pe);
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarTopCRTHitPe([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=30&&hit.plane<=39){
        rets.push_back(hit.pe);
      }
    }
    return rets;
  });
  // - Pos
  const SpillMultiVar spillvarEastWestCRTHitPosX([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=40&&hit.plane<=45){
        rets.push_back(hit.position.x);
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarFlashPosX([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      rets.push_back(opflash.center.x);
    }
    return rets;
  });
  // - OpFlash
  const SpillMultiVar spillvarOpFlashPeakToFirstTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      rets.push_back( 1000.*(opflash.time - opflash.firsttime) );
    }
    return rets;
  });
  // - CRTHit
  const SpillMultiVar spillvarTopCRTHitTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=30&&hit.plane<=39){
        double this_crttime = sr->hdr.ismc ? hit.t0 : hit.t1;
        rets.push_back(this_crttime);
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarSideCRTHitTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=40&&hit.plane<=49){
        double this_crttime = sr->hdr.ismc ? hit.t0 : hit.t1;
        rets.push_back(this_crttime);
      }
    }
    return rets;
  });
  // - PMT-CRT matching
  const SpillMultiVar spillvarTopCRTPMTTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = ICARUSCRTPMTMatching::spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : intimeTimes){
      std::vector<int> crtHitIdices = ICARUSCRTPMTMatching::cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
      for(const auto& crtHitIdx:crtHitIdices){
        double this_crttime = sr->hdr.ismc ? sr->crt_hits.at(crtHitIdx).t0 : sr->crt_hits.at(crtHitIdx).t1;
        rets.push_back( this_crttime - opt );
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarSideCRTPMTTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = ICARUSCRTPMTMatching::spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : intimeTimes){
      std::vector<int> crtHitIdices = ICARUSCRTPMTMatching::cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 1);
      for(const auto& crtHitIdx:crtHitIdices){
        double this_crttime = sr->hdr.ismc ? sr->crt_hits.at(crtHitIdx).t0 : sr->crt_hits.at(crtHitIdx).t1;
        rets.push_back( this_crttime - opt );
      }
    }
    return rets;
  });

  // - WW test
  const SpillMultiVar spillvarWWTPCTrackEndX([](const caf::SRSpillProxy *sr){
    vector<double> rets;
    for(const auto& slc: sr->slc){
      for(const auto& pfp: slc.reco.pfp){
        const auto& trk = pfp.trk;
        if( !isnan(trk.start.x) && !isnan(trk.start.y) && !isnan(trk.start.z)
            && !isnan(trk.end.x) && !isnan(trk.end.x) && !isnan(trk.end.x) ){
          rets.push_back( trk.start.x );
          rets.push_back( trk.end.x );
        }
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarWWTPCTrackEndY([](const caf::SRSpillProxy *sr){
    vector<double> rets;
    for(const auto& slc: sr->slc){
      for(const auto& pfp: slc.reco.pfp){
        const auto& trk = pfp.trk;
        if( !isnan(trk.start.x) && !isnan(trk.start.y) && !isnan(trk.start.z)
            && !isnan(trk.end.x) && !isnan(trk.end.x) && !isnan(trk.end.x) ){
          if( trk.start.x>220. && trk.start.x<340. ){
            rets.push_back( trk.start.y );
            rets.push_back( trk.end.y );
          }
        }
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarWWTPCTrackEndZ([](const caf::SRSpillProxy *sr){
    vector<double> rets;
    for(const auto& slc: sr->slc){
      for(const auto& pfp: slc.reco.pfp){
        const auto& trk = pfp.trk;
        if( !isnan(trk.start.x) && !isnan(trk.start.y) && !isnan(trk.start.z)
            && !isnan(trk.end.x) && !isnan(trk.end.x) && !isnan(trk.end.x) ){
          if( trk.start.x>220. && trk.start.x<340. ){
            rets.push_back( trk.start.z );
            rets.push_back( trk.end.z );
          }
        }
      }
    }
    return rets;
  });

  // Var
  // - GENIE interaction code
  // - https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360
  const Var varGENIEIntCode([](const caf::SRSliceProxy* slc) -> int {
    if(slc->truth.index >= 0.f){
      if(slc->truth.genie_mode<0) return -1;
      else if(slc->truth.genie_mode>13) return 14;
      else return slc->truth.genie_mode;
    }
    else{
      return -2;
    }
  });
  // - Truth interaction
  const Var varNeutrinoTruthE([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.E)) return -999.;
    else return slc->truth.E;
  });
  const Var varNuDirectionX([](const caf::SRSliceProxy* slc) ->int {
      double this_x = slc->truth.prod_vtx.x/100.;
      double this_y = slc->truth.prod_vtx.y/100.;
      double this_z = slc->truth.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z).Unit();
      this_coord *= -1.;
      return this_coord.X();
  });
  const Var varNuDirectionY([](const caf::SRSliceProxy* slc) ->int {
      double this_x = slc->truth.prod_vtx.x/100.;
      double this_y = slc->truth.prod_vtx.y/100.;
      double this_z = slc->truth.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z).Unit();
      this_coord *= -1.;
      return this_coord.Y();
  });
  const Var varNuDirectionZ([](const caf::SRSliceProxy* slc) ->int {
      double this_x = slc->truth.prod_vtx.x/100.;
      double this_y = slc->truth.prod_vtx.y/100.;
      double this_z = slc->truth.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z).Unit();
      this_coord *= -1.;
      return this_coord.Z();
  });
  const Var varTruthQ2([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.Q2)) return -999.;
    else return slc->truth.Q2;
  });
  const Var varTruthq0_lab([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.q0_lab)) return -999.;
    else return slc->truth.q0_lab;
  });
  const Var varTruthmodq_lab([](const caf::SRSliceProxy* slc) -> double {
    double Q2 = varTruthQ2(slc);
    double q0_lab = varTruthq0_lab(slc);
    if(isnan(Q2)||isnan(q0_lab)) return -999.;
    else{
      return sqrt(Q2*Q2+q0_lab*q0_lab);
    }
  });
  const Var varTruthW([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.w)) return -999.;
    else return slc->truth.w;
  });
  // - Slice var
  const Var varCountSlice([](const caf::SRSliceProxy* slc) ->int {
    return 0.;
  });
  const Var varIsClearCosmic([](const caf::SRSliceProxy* slc) ->int {
    if(slc->is_clear_cosmic) return 1;
    else return 0;
  });
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
  // - NuID
  const Var varCRLongestTrackDirY([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.crlongtrkdiry;
  });
  // - Long-enough tracks
  const MultiVar varLongTrackDirectionY([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> rets;
    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){

      bool IsTrack = slc->reco.pfp.at(i).trackScore > 0.5;
      if(!IsTrack) continue;
      const auto& trk = slc->reco.pfp.at(i).trk;
      if(isnan(trk.len)) continue;
      if(trk.len>50.) rets.push_back(trk.dir.y);

    }
    return rets;
  });
  // - Primary tracks
  const MultiVar PrimaryTrackIndices([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> rets;
    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
      const auto& pfp = slc->reco.pfp.at(i_pfp);
      if(pfp.trackScore<0.5) continue;
      const auto& trk = pfp.trk;
      if(isnan(trk.start.x)) continue;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      if (Atslc < 10. && pfp.parent_is_primary){
        rets.push_back(i_pfp);
      }
    }
    return rets;
  });
  const Var NPrimaryTracks([](const caf::SRSliceProxy* slc) -> double {
    return PrimaryTrackIndices(slc).size();
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

      if(isnan(trk.len)) continue;
      if(isnan(trk.end.x)) continue;

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
