#include "sbnana/SBNAna/CRTPMTMatching/ICARUSCRTPMTMatching_Vars.h"

using namespace std;
using namespace ana;

namespace ICARUSCRTPMTMatching{

  // - PMT-CRT matching
  const SpillMultiVar spillvarOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      double this_oft = opflash.firsttime;
      //if(sr->hdr.ismc) this_oft += -55.1e-3; // https://github.com/SBNSoftware/icaruscode/blob/v09_64_01/icaruscode/CRT/CrtOpHitMatchAnalysis_module.cc#L187
      if(sr->hdr.ismc) this_oft += -43.0e-3; // https://github.com/SBNSoftware/icaruscode/blob/v09_64_01/icaruscode/CRT/CrtOpHitMatchAnalysis.fcl#L60
      rets.push_back(this_oft);
    }
    return rets;
  });
  const SpillMultiVar spillvarValidOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      double this_oft = opflash.firsttime;
      //if(sr->hdr.ismc) this_oft += -55.1e-3; // https://github.com/SBNSoftware/icaruscode/blob/v09_64_01/icaruscode/CRT/CrtOpHitMatchAnalysis_module.cc#L187
      if(sr->hdr.ismc) this_oft += -43.0e-3; // https://github.com/SBNSoftware/icaruscode/blob/v09_64_01/icaruscode/CRT/CrtOpHitMatchAnalysis.fcl#L60
      if(opflash.onbeamtime){
        rets.push_back(this_oft); // TODO I'm using a hacked version of onbeamtime.. in my pricate samples, onbeamtime = "has >=5 OpHits with ADC>400"
      }
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
      std::vector<int> crtHitIdices = cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
      for(const auto& crtHitIdx:crtHitIdices){
        double this_crttime = sr->hdr.ismc ? sr->crt_hits.at(crtHitIdx).t0 : sr->crt_hits.at(crtHitIdx).t1;
        rets.push_back( this_crttime - opt );
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

  const SpillMultiVar spillvarMatchID2_CRTHitPosXs([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    if(intimeTimes.size()>0){
      for(const auto& opt : intimeTimes){
        int MatchID = cpmt.GetMatchID(opt, sr->crt_hits);
        if(MatchID==2){
          vector<int> crtHitIdices = cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
          for(const auto& crtHitIdx:crtHitIdices){
            const auto& crthit_matched = sr->crt_hits.at(crtHitIdx);
            rets.push_back( crthit_matched.position.x );
          }
        }
      }
    }
    return rets;
  });



} // end namespace ICARUSCRTPMTMatching
