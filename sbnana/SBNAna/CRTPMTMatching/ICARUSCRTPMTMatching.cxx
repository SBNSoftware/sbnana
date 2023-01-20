#include "sbnana/SBNAna/CRTPMTMatching/ICARUSCRTPMTMatching.h"

using namespace ana;

namespace ICARUSCRTPMTMatching{

  CRTPMTMatchingTool::CRTPMTMatchingTool(){
    std::cout << "[CRTPMTMatchingTool::CRTPMTMatchingTool] called" << std::endl;
    UseTS0 = false;
  }

  CRTPMTMatchingTool& CRTPMTMatchingTool::Instance(){
    static CRTPMTMatchingTool cpmt;
    return cpmt;
  }

  void CRTPMTMatchingTool::SetGateType(GateType gt) const {
    GT = gt;
    if(abs(GT)==1){
      timecut_min = 0.;
      timecut_max = 2.2;
    }
    else if(abs(GT)==2){
      timecut_min = 0.;
      timecut_max = 10.1;
    }
    else{
      std::cout << "[CRTPMTMatchingTool] Wrong gate type = " << gt << std::endl;
      abort();
    }

    std::cout << "[CRTPMTMatchingTool::SetGateType] GT = " << GT << ", timecut_min = " << timecut_min << ", timecut_max = " << timecut_max << std::endl;

  }

  void CRTPMTMatchingTool::SetInTimeRange(double t_min, double t_max) const {
    timecut_min = t_min;
    timecut_max = t_max;
    std::cout << "[CRTPMTMatchingTool::SetInTimeRange] timecut_min = " << timecut_min << ", timecut_max = " << timecut_max << std::endl;
  }

  bool CRTPMTMatchingTool::IsInTime(double t_gate) const{
    //std::cout << "timecut_min = " << timecut_min << ", t_gate = " << t_gate << ", timecut_max = " << timecut_max << std::endl;
    return ( timecut_min<=t_gate && t_gate<=timecut_max );
  }

  int CRTPMTMatchingTool::GetMatchedCRTHitIndex(
    double opt,
    const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
    int mode) const{

    double mindiff = std::numeric_limits<double>::max();
    int ret=-1;
    for(size_t i=0; i<crt_hits.size(); i++){
      const auto& hit = crt_hits.at(i);
      if(hit.plane>=30 && hit.plane<=34){
        double crtt = UseTS0 ? hit.t0 : hit.t1;
        double this_diff = crtt-opt;
        if(fabs(this_diff)<fabs(mindiff)){
          mindiff = this_diff;
          ret = i;
        }
      }
    }

    return ret;

  }

  int CRTPMTMatchingTool::GetMatchID(
    double opt,
    const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits) const{

    int hasCRTHit = 0;
    static double interval = 0.1;
    int topen = 0, topex = 0, sideen = 0, sideex = 0;
    for(size_t i=0; i<crt_hits.size(); i++){
      const auto& hit = crt_hits.at(i);
      double crtt = UseTS0 ? hit.t0 : hit.t1;
      double tof = crtt - opt;

      if(tof<0 && abs(tof)<interval){
        if(hit.plane > 36){
          sideen++;
        }
        else{
          topen++;
        }
      }
      else if(tof>=0 && abs(tof)<interval){
        if(hit.plane > 36){
          sideex++;
        }
        else{
          topex++;
        }
      }

    }

    // hasCRTHit = 0, no matched CRT
    // hasCRTHit = 1, 1 entering from Top CRT
    // hasCRTHit = 2, 1 entering from Side CRT
    // hasCRTHit = 3, 1 entering from Top and exiting to Side CRT
    // hasCRTHit = 4, No entering; 1 exiting to top
    // hasCRTHit = 5, No entering, 1 exiting to side
    // hasCRTHit = 6, Multiple entering
    // hasCRTHit = 7, Multiple entering and exiting to side
    // hasCRTHit = 8, all other cases

    if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0)
      hasCRTHit = 0;
    else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0)
      hasCRTHit = 1;
    else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0)
      hasCRTHit = 2;
    else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1)
      hasCRTHit = 3;
    else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0)
      hasCRTHit = 4;
    else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1)
      hasCRTHit = 5;
    else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0)
      hasCRTHit = 6;
    else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1)
      hasCRTHit = 7;
    else
      hasCRTHit = 8;

    return hasCRTHit;

  }

  bool CRTPMTMatchingTool::IsNegativeTOF(double timediff) const{
    return (timediff>-0.1 && timediff<0.);
  }

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
  const SpillMultiVar spillvarCRTPMTTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : intimeTimes){
      int crtHitIdx = cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
      if(crtHitIdx>=0){
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

  // - CRTPMT matching

  const SpillCut spillcutHasValidFlash([](const caf::SRSpillProxy* sr){
    vector<double> validTimes = spillvarValidOpFlashTime(sr);
    return validTimes.size()>0;
  });
  const SpillCut spillcutHasInTimeFlash([](const caf::SRSpillProxy* sr){
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr); // Valid and InTime
    return intimeTimes.size()>0;
  });
  const SpillCut spillcutCRTPMTCosmicByID([](const caf::SRSpillProxy* sr){

    // A spill is considered as "Cosmic" when all intime flahses are entering
    // If a spill has anything not entering, that shuold not be considered as a cosmic

    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    if(intimeTimes.size()>0){

      bool IsAllEnteringInThisSpill = true;
      for(const auto& opt : intimeTimes){
        int crtHitIdx = cpmt.GetMatchID(opt, sr->crt_hits);
        bool IsThisFlashEntering = (crtHitIdx==1 || crtHitIdx==2 || crtHitIdx==3 || crtHitIdx==6 || crtHitIdx==7);
        if(!IsThisFlashEntering) IsAllEnteringInThisSpill = false;
      }
      return IsAllEnteringInThisSpill;

    }
    else{
      return false;
    }
  });

  const SpillCut spillcutCRTPMTHasNegativeTOF([](const caf::SRSpillProxy* sr){
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    if(intimeTimes.size()>0){
      for(const auto& opt : intimeTimes){
        int crtHitIdx = cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
        if(crtHitIdx>=0){
          bool IsNegativeTOF = cpmt.IsNegativeTOF( sr->crt_hits.at(crtHitIdx).t1 - opt );
          if(IsNegativeTOF) return true;
        }
      }
      return false;
    }
    else{
      return false;
    }
  });
  const SpillCut spillcutCRTPMTAllNegativeTOF([](const caf::SRSpillProxy* sr){
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    if(intimeTimes.size()>0){
      for(const auto& opt : intimeTimes){
        int crtHitIdx = cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
        if(crtHitIdx>=0){
          bool IsNegativeTOF = cpmt.IsNegativeTOF( sr->crt_hits.at(crtHitIdx).t1 - opt );
          if(!IsNegativeTOF) return false;
        }
      }
      return true;
    }
    else{
      return false;
    }
  });


} // END namespace
