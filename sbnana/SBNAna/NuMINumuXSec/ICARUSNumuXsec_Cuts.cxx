#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  // SpillCut

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

  // Cut (slice)

  // - FV

  const Cut cutRecoVertexFV([](const caf::SRSliceProxy* slc) {

    if( !isnan(slc->vertex.x) ) return fv.isContained(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    else return false;

  });

  // - Flash matching

  const Cut cutFMScore([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 12.0 && slc->fmatch.score >= 0 );
  });

  const Cut cutFMTime([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.time) && slc->fmatch.time>=-0.2 && slc->fmatch.time<=9.9 );
  });

}
