#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  // SpillCut

  // - CRTPMT matching

  const SpillCut spillcutHasValidFlash([](const caf::SRSpillProxy* sr){
    vector<double> validTimes = ICARUSCRTPMTMatching::spillvarValidOpFlashTime(sr);
    return validTimes.size()>0;
  });
  const SpillCut spillcutHasInTimeFlash([](const caf::SRSpillProxy* sr){
    vector<double> intimeTimes = ICARUSCRTPMTMatching::spillvarInTimeOpFlashTime(sr); // Valid and InTime
    return intimeTimes.size()>0;
  });
  const SpillCut spillcutHasOneTimeFlash([](const caf::SRSpillProxy* sr){
    vector<double> intimeTimes = ICARUSCRTPMTMatching::spillvarInTimeOpFlashTime(sr); // Valid and InTime
    return intimeTimes.size()==1;
  });

  const SpillCut spillcutAllEnteringByCRTPMT([](const caf::SRSpillProxy* sr){

    // A spill is considered as "Cosmic" when all intime flahses are entering
    // If a spill has anything not entering, that shuold not be considered as a cosmic

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

    vector<double> intimeOpFlashIDs = ICARUSCRTPMTMatching::spillvarCRTPMTMatchingID(sr);
    if(intimeOpFlashIDs.size()>0){

      bool IsAllEnteringInThisSpill = true;
      for(const auto& opID : intimeOpFlashIDs){
        //bool IsThisFlashEntering = (opID==1 || opID==2 || opID==3 || opID==4 || opID==5 || opID==6 || opID==7 || opID==8 || opID==9 || opID==10 || opID==11);
        bool IsThisFlashEntering = (opID>=1 && opID<=11);
        if(!IsThisFlashEntering) IsAllEnteringInThisSpill = false;
      }
      return IsAllEnteringInThisSpill;

    }
    else{
      return false;
    }

  });

  const SpillCut spillcutHasCRTPMTDirtByID([](const caf::SRSpillProxy* sr){

    // A spill is considered as "Cosmic" when all intime flahses are entering
    // If a spill has anything not entering, that shuold not be considered as a cosmic

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

    vector<double> intimeOpFlashIDs = ICARUSCRTPMTMatching::spillvarCRTPMTMatchingID(sr);
    if(intimeOpFlashIDs.size()>0){

      for(const auto& opID : intimeOpFlashIDs){
        bool IsThisFlashDirt = (opID==3 || opID==4 || opID==5 || opID==6);
        if(IsThisFlashDirt) return true;
      }
      return false;

    }
    else{
      return false;
    }

  });

  // - Spill nu
  const SpillCut spillcutHasTrueMuContained([](const caf::SRSpillProxy* sr){

    for(const auto& _nu : sr->mc.nu){
      for(const auto& _prim: _nu.prim){
        if( abs(_prim.pdg)==13 && _prim.contained ) return true;
      }
    }
    return false;

  });
  const SpillCut spillcutHasTrueMuExiting([](const caf::SRSpillProxy* sr){
    
    for(const auto& _nu : sr->mc.nu){
      for(const auto& _prim: _nu.prim){
        if( abs(_prim.pdg)==13 && !(_prim.contained) ) return true;
      }
    }
    return false;

  });

  // Cut (slice)

  // - Interaction
  //   - Neutrino flavor
  const Cut cutIsNuMu([](const caf::SRSliceProxy* slc) {
    return (kIsNuSlice(slc) && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ));
  });
  const Cut cutIsNuE([](const caf::SRSliceProxy* slc) {
    return (kIsNuSlice(slc) && ( slc->truth.pdg == 12 || slc->truth.pdg == -12 ));
  });
  //   - CC vs NC
  const Cut cutIsCC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.iscc );
  });
  const Cut cutIsNC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.isnc );
  });
  //   - GENIE Interaction code
  const Cut cutIsQE([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==0;
  });
  const Cut cutIsRes([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==1;
  });
  const Cut cutIsDIS([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==2;
  });
  const Cut cutIsCoh([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==3;
  });
  const Cut cutIsCohElastic([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==4;
  });
  const Cut cutIsElectronScattering([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==5;
  });
  const Cut cutIsIMDAnnihilation([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==6;
  });
  const Cut cutIsInverseBetaDecay([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==7;
  });
  const Cut cutIsGlashowResonance([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==8;
  });
  const Cut cutIsAMNuGamma([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==9;
  });
  const Cut cutIsMEC([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==10;
  });
  const Cut cutIsDiffractive([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==11;
  });
  const Cut cutIsEM([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==12;
  });
  const Cut cutIsWeakMix([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==13;
  });
  const Cut cutIsUnknownInteractionType1([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==-1;
  });
  const Cut cutIsUnknownInteractionType2([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)>13;
  });
  const Cut cutIsUnknownInteractionType3([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)<-1;
  });
  //   - NuMu-CC categories
  const Cut cutIsNuMuCC([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMu(slc) && cutIsCC(slc) );
  });
  const Cut cutIsNuMuCCQE([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsQE(slc) );
  });
  const Cut cutIsNuMuCCRes([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsRes(slc) );
  });
  const Cut cutIsNuMuCCMEC([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsMEC(slc) );
  });
  const Cut cutIsNuMuCCDIS([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsDIS(slc) );
  });
  const Cut cutIsNuMuCCCoh([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsCoh(slc) );
  });
  const Cut cutIsNuMuCCCohElastic([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsCohElastic(slc) );
  });
  //   - NuMu-NC categories
  const Cut cutIsNuMuNC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.isnc && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ) );
  });
  //   - NuE
  const Cut cutIsNuECC([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuE(slc) && cutIsCC(slc) );
  });

  // - FV

  const Cut cutRecoVertexFV([](const caf::SRSliceProxy* slc) {

    if( !isnan(slc->vertex.x) ) return fv.isContained(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    else return false;

  });
  const Cut cutTruthVertexActiveVolume([](const caf::SRSliceProxy* slc) {
    if( !isnan(slc->truth.position.x) ) return av.isContained(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
    else return false;
  });
  const Cut cutRFiducial([](const caf::SRSliceProxy* slc) {

    if( !isnan(slc->vertex.x) ) return fv.isContained(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    else return false;

  });

  // - Pandora based ClearCosmic

  const Cut cutNotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !(slc->is_clear_cosmic);
  });
  const Cut cutCRLongestTrackDirY([](const caf::SRSliceProxy* slc) {
    return (slc->nuid.crlongtrkdiry>-0.7);
  });
  const Cut cutCRLongestTrackDirYHard([](const caf::SRSliceProxy* slc) {
    return (slc->nuid.crlongtrkdiry>-0.4);
  });

  // - Flash matching

  const Cut cutFMScore([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 12.0 && slc->fmatch.score >= 0 );
  });

  const Cut cutFMTime([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.time) && slc->fmatch.time>=-0.2 && slc->fmatch.time<=9.9 );
  });

}
