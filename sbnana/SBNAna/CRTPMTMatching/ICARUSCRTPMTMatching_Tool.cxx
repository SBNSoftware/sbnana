#include "sbnana/SBNAna/CRTPMTMatching/ICARUSCRTPMTMatching_Tool.h"

using namespace ana;

namespace ICARUSCRTPMTMatching{

CRTPMTMatchingTool::CRTPMTMatchingTool(){
  std::cout << "[CRTPMTMatchingTool::CRTPMTMatchingTool] called" << std::endl;
  UseTS0 = false;
  Debug = false;
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

// mode = 0 : Top+Side
// mode = 1 : Top only
// mode = 2 : Side only
std::vector<int> CRTPMTMatchingTool::GetMatchedCRTHitIndex(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
  int mode) const{

  static double interval = 0.1;

  //double mindiff = std::numeric_limits<double>::max();
  std::vector<int> rets = {};
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& hit = crt_hits.at(i);
    bool IsTopCRT = hit.plane>=30 && hit.plane<=39;
    bool IsSideCRT = hit.plane>=40 && hit.plane<=49;
    bool crtSelection = false;
    if(mode==0) crtSelection = IsTopCRT||IsSideCRT;
    if(mode==1) crtSelection = IsTopCRT;
    if(mode==2) crtSelection = IsSideCRT;
    if( crtSelection ){
      double crtt = UseTS0 ? hit.t0 : hit.t1;
      double this_diff = crtt-opt;
      if(abs(this_diff)<interval) rets.push_back(i);
/*
      if(fabs(this_diff)<fabs(mindiff)){
        mindiff = this_diff;
        ret = i;
      }
*/
    }
  }

  return rets;

}

int CRTPMTMatchingTool::GetMatchID(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
  int mode) const{

  int hasCRTHit = -1;
  static double interval = 0.1;
  int TopEn = 0, TopEx = 0, SideEn = 0, SideEx = 0;
  int SideSouthEn = 0;
  int SideEastEn = 0;
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& hit = crt_hits.at(i);
    double crtt = UseTS0 ? hit.t0 : hit.t1;
    double tof = crtt - opt;

    if(tof<0 && abs(tof)<interval){
      if(hit.plane >= 40 && hit.plane<=49){
        SideEn++;
        if(hit.plane==46){
          SideSouthEn++;
        }
        if(hit.plane>=43&&hit.plane<=45){
          SideEastEn++;
        }
      }
      else if(hit.plane >= 30 && hit.plane <= 39){
        TopEn++;
      }
      else{
      }
    }
    else if(tof>=0 && abs(tof)<interval){
      if(hit.plane >= 40 && hit.plane<=49){
        SideEx++;
      }
      else if(hit.plane >= 30 && hit.plane <= 39){
        TopEx++;
      }
      else{
      }
    }

  }

  int OtherSideEn = SideEn - SideSouthEn - SideEastEn;

  // Top+Side
  if(mode==0){

    // -- No matching
    // hasCRTHit = -1, no matched CRT
    // -- 0 entering, but somthing exiting cases (e.g., exiting muon from nu)
    // hasCRTHit = 0, 0 entering, something exiting
    // -- Top-entering cases (Cosmic)
    // hasCRTHit = 1, 1 entering from Top CRT, nothing else
    // hasCRTHit = 2, 1 entering from Top CRT, >=1 exiting to side
    // -- South-entering cases (BNB or NuMI DIRT)
    // hasCRTHit = 3, 1 entering from South-Side CRT (no East-Side entering), nothing else
    // hasCRTHit = 4, 1 entering from South-Side CRT (no East-Side entering), exiting to somewhere
    // -- East-entering cases (NuMI DIRT)
    // hasCRTHit = 5, 1 entering from East-Side CRT (no South-Side entering), nothing else
    // hasCRTHit = 6, 1 entering from Esat-Side CRT (no South-Side entering), exiting to somewhere
    // -- Other side wall entering cases (Cosmic)
    // hasCRTHit = 7, 1 entering from Other-Side CRT, nothing else
    // hasCRTHit = 8, 1 entering from Other-Side CRT, exiting to somewhere
    // -- Multiple top-entering
    // hasCRTHit = 9, >=1 entering from Top CRT, nothing else
    // hasCRTHit = 10, >=1 entering from Top CRT, >=1 exiting to side
    // --- Top-Side entering
    // hasCRTHit = 11, >=1 entering from Top CRT, >=1 entering from side
    // --- TEST
    // hasCRTHit = 12, >=1 entering top and >=1 exiting top

    if(TopEn==0 && SideEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = -1;

    else if(TopEn==0 && SideEn==0 && TopEx+SideEx>=1)
      hasCRTHit = 0;

    else if(TopEn==1 && SideEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 1;
    else if(TopEn==1 && SideEn==0 && TopEx==0 && SideEx>=1)
      hasCRTHit = 2;

    else if(TopEn==0 && SideSouthEn==1 && SideEastEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 3;
    else if(TopEn==0 && SideSouthEn==1 && SideEastEn==0 && TopEx+SideEx>=1)
      hasCRTHit = 4;

    else if(TopEn==0 && SideEastEn==1 && SideSouthEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 5;
    else if(TopEn==0 && SideEastEn==1 && SideSouthEn==0 && TopEx+SideEx>=1)
      hasCRTHit = 6;

    else if(TopEn==0 && OtherSideEn==1 && TopEx==0 && SideEx==0)
      hasCRTHit = 7;
    else if(TopEn==0 && OtherSideEn==1 && TopEx+SideEx>=1)
      hasCRTHit = 8;

    else if(TopEn>=1 && SideEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 9;
    else if(TopEn>=1 && SideEn==0 && TopEx==0 && SideEx>=1)
      hasCRTHit = 10;

    else if(TopEn>=1 && SideEn>=1)
      hasCRTHit = 11;

    else if(TopEn>=1 && TopEx>=1)
      hasCRTHit = 12;

    else{
      if(Debug) printf("(TopEn, TopEx, SideEn, SideEx) = (%d, %d, %d, %d)\n",TopEn, TopEx, SideEn, SideEx);
      hasCRTHit = 13;
    }

  }
  // Top only
  else if(mode==1){

    // hasCRTHit = 0, no matched CRT
    // hasCRTHit = 1, One top entering (no Top-exiting)
    // hasCRTHit = 2, One top exiting (no Top-entering)
    // hasCRTHit = 3, Multiple top entering (no top-exiting)
    // hasCRTHit = 4, Multiple top exiting (no top-entering)
    // hasCRTHit = 5, All other cases; entering>=1 && exiting>=1

    if(TopEn==0 && TopEx==0)
      hasCRTHit = 0;
    else if(TopEn==1 && TopEx==0)
      hasCRTHit = 1;
    else if(TopEn==0 && TopEx==1)
      hasCRTHit = 2;
    else if(TopEn>=1 && TopEx==0)
      hasCRTHit = 3;
    else if(TopEn==0 && TopEx>=1)
      hasCRTHit = 4;
    else
      hasCRTHit = 5;

  }
  // side only
  else if(mode==2){

    // hasCRTHit = 0, no matched CRT
    // hasCRTHit = 1, One side entering (no side-exiting)
    // hasCRTHit = 2, One side exiting (no side-entering)
    // hasCRTHit = 3, Multiple side entering (no side-exiting)
    // hasCRTHit = 4, Multiple side exiting (no side-entering)
    // hasCRTHit = 5, All other cases; entering>=1 && exiting>=1

    if(SideEn==0 && SideEx==0)
      hasCRTHit = 0;
    else if(SideEn==1 && SideEx==0)
      hasCRTHit = 1;
    else if(SideEn==0 && SideEx==1)
      hasCRTHit = 2;
    else if(SideEn>=1 && SideEx==0)
      hasCRTHit = 3;
    else if(SideEn==0 && SideEx>=1)
      hasCRTHit = 4;
    else
      hasCRTHit = 5;

  }
  else{

  }

  return hasCRTHit;

}

bool CRTPMTMatchingTool::IsNegativeTOF(double timediff) const{
  return (timediff>-0.1 && timediff<0.);
}

} // END namespace
