#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/SBNAna/Cuts/BNBQualityCuts_Jan2025.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <iostream>

namespace ana { // use namespace ana for entire file

/*
   Here are the hard-coded physically meaningful bounds for BNB quality
   SpillCuts. 
*/

/*
   While there is sometimes discernable structure in a given BNB quality 
   variable (e.g. TOR860 and TOR875 stay around 4.5E12 POT when only the BNB is
   on and 3.0E12 POT when both the BNB and NuMI beams are on), exact behavior of
   such variables is not known in detail as a function of time. Thus, the latest 
   (read 'final') version of a given cut will be the most stringent cut we can 
   make with a limited knowledge of the detailed behavior of that cut.
*/

const double TOR860_LB = +100e9; // units: POT
const double TOR875_LB = +100e9; 

const double LM875A_LB = +1e-2; // units: rad/s
const double LM875B_LB = +1e-2;
const double LM875C_LB = +1e-2;

// SMITHJA: I will have to get in touch with (Tom from) the External Beams group
//          at Fermilab to get final cuts on these beam position monitor (BPM)
//          variables since they change from run to run. If realtime values for 
//          those optimal BPM variables are not accessible (which is the likely
//          case), the BPM cuts will likely be hard-coded here. 
const double HP875_LB = -0.25;  const double HP875_UB = +0.25; // units: mm
const double VP875_LB = -0.20;  const double VP875_UB = +0.30;
const double HPTG1_LB = -0.25;  const double HPTG1_UB = +0.25;
const double VPTG1_LB = -0.20;  const double VPTG1_UB = +0.30;
const double HPTG2_LB = -0.25;  const double HPTG2_UB = +0.25;
const double VPTG2_LB = -0.20;  const double VPTG2_UB = +0.30;

/*
   The BTJT2 variable is somewhat unreliable to use as a spill-by-spill cut.
   BTJT2 measures the temperature of the BNB near the target. Since runs can be
   started/stopped more quickly than the thermal fluctuations of BTJT2, there
   can be some amount of "cross contamination" between spills. For this reason
   BTJT2 is deemed to have diminishing returns compared to other BNB quality
   variables listed. You may employ cuts on BTJT2 at your own peril.

const double BTJT2_LB = 999;  const double BTJT2_UB = 999; // units: deg C
*/

const double THCURR_LB = +173;  const double THCURR_UB = +177; // units: kA

/* SMITHJA:
   Further cuts on THCURR could/should be made since it was discovered that
   the horn current (i.e. THCURR) is not always centered around 174 kA. Further
   work needed to understand horn current for different runs. -- Jacob Smith

const double THCURR_LB = +173;  const double THCURR_UB = +175; // units: kA
*/
//-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



struct SpillInfo { // contains information about a spill of the BNB
  // These are the variables that measure the quality of a spill in the BNB.
  float TOR860, TOR875, \
        LM875A, LM875B, LM875C, \
        HP875, VP875, \
        HPTG1, VPTG1, \
        HPTG2, VPTG2, \
//        BTJT2,
        THCURR;
  
  // These are used to identify/match spills in the global vector of 
  // SpillInfo objects (see GLOBAL_SPILL_INFO below).
  unsigned int run, event;
};

// vector containing all of the information for spills of the BNB:
static std::vector<SpillInfo> GLOBAL_SPILL_INFO; 


/*
   This function sets all of the values for a SpillInfo object (see above).

Params:
  * run: integer identifier for the run/event
  * info: this contains all of the BNB monitoring info for a given spill 
          (including the event number, which is also used to identify a
           run/event)
*/
void fill_spill_info( unsigned int run, const caf::Proxy<caf::SRBNBInfo> &info){
  GLOBAL_SPILL_INFO.emplace_back();

  GLOBAL_SPILL_INFO.back().run = run;
  GLOBAL_SPILL_INFO.back().event = info.event;

  GLOBAL_SPILL_INFO.back().TOR860 = info.TOR860;
  GLOBAL_SPILL_INFO.back().TOR875 = info.TOR875;
  GLOBAL_SPILL_INFO.back().LM875A = info.LM875A;
  GLOBAL_SPILL_INFO.back().LM875B = info.LM875B;
  GLOBAL_SPILL_INFO.back().LM875C = info.LM875C;
  GLOBAL_SPILL_INFO.back().HP875  = info.HP875;
  GLOBAL_SPILL_INFO.back().VP875  = info.VP875;
  GLOBAL_SPILL_INFO.back().HPTG1  = info.HPTG1;
  GLOBAL_SPILL_INFO.back().VPTG1  = info.VPTG1;
  GLOBAL_SPILL_INFO.back().HPTG2  = info.HPTG2;
  GLOBAL_SPILL_INFO.back().VPTG2  = info.VPTG2;
//  GLOBAL_SPILL_INFO.back().BTJT2  = info.BTJT2;
  GLOBAL_SPILL_INFO.back().THCURR = info.THCURR;
}

const SpillCut kTOR860([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){                           // non-zero nbnb indicates new subrun...
        GLOBAL_SPILL_INFO.clear();                   // ...so we get rid of spills tied to last subrun...
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ // ...and fill in info for this subrun's spills
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; // ensure event matching
        // make preliminary cut on variable of interest --v
        if( match.TOR860 >= TOR860_LB) return true; 
    }
    return false; // a spill does not pass a cut if 
                  // there is no event-matched data
});

/*
   NOTE: For brevity of commenting, the preliminary SpillCuts below are not 
         commented; however, they are the same as kTOR860 (see above) but with
         different variables of interest, e.g. TOR875, LM875B, THCURR, etc.
*/

const SpillCut kTOR875([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){                               
        GLOBAL_SPILL_INFO.clear();                   
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.TOR875 >= TOR875_LB) return true; 
    }
    return false; 
});

const SpillCut kLM875A([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){                               
        GLOBAL_SPILL_INFO.clear();                   
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.LM875A >= LM875A_LB) return true; 
    }
    return false; 
});

const SpillCut kLM875B([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){                               
        GLOBAL_SPILL_INFO.clear();                   
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.LM875B >= LM875B_LB) return true; 
    }
    return false; 
});

const SpillCut kLM875C([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){                               
        GLOBAL_SPILL_INFO.clear();                   
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.LM875C >= LM875C_LB) return true; 
    }
    return false; 
});

const SpillCut kHP875([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.HP875 >= HP875_LB && match.HP875 <= HP875_UB) return true;
    }
    return false; 
});

const SpillCut kVP875([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.VP875 >= VP875_LB && match.VP875 <= VP875_UB) return true;
    }
    return false; 
});

const SpillCut kHPTG1([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.HPTG1 >= HPTG1_LB && match.HPTG1 <= HPTG1_UB) return true;
    }
    return false; 
});

const SpillCut kVPTG1([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.VPTG1 >= VPTG1_LB && match.VPTG1 <= VPTG1_UB) return true;
    }
    return false; 
});

const SpillCut kHPTG2([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.HPTG2 >= HPTG2_LB && match.HPTG2 <= HPTG2_UB) return true;
    }
    return false; 
});

const SpillCut kVPTG2([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.VPTG2 >= VPTG2_LB && match.VPTG2 <= VPTG2_UB) return true;
    }
    return false; 
});
/*
const SpillCut kBTJT2([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.BTJT2 >= BTJT2_LB && match.BTJT2 < BTJT2_UB) return true;
    }
    return false; 
});
*/
const SpillCut kTHCURR([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.THCURR >= THCURR_LB && match.THCURR <= THCURR_UB) return true;
    }
    return false; 
});

/* -------------------------------------- \
//                                        |
// ^--- INDIVIDUAL CUTS ; MASTER CUT ---v |
//                                        |
// ------------------------------------- */

const SpillCut kBNBQuality_noBPMs_10Feb2025([](const caf::SRSpillProxy* sr) {
    if( sr.hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hrd.bnbinfo){
            fill_spill_info( sr->hdr.run, bnbinfo);
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue;
        if( match.TOR860

        if( match.event != sr->hdr.evt) continue; 
        if( match.TOR860 >= TOR860_LB && \
              match.TOR875 >= TOR875_LB && \
              match.LM875A >= LM875A_LB && \
              match.LM875B >= LM875B_LB && \
              match.LM875C >= LM875C_LB && \
//              match.HP875  >= HP875_LB  && match.HP875 <= HP875_UB  && 
//              match.VP875  >= VP875_LB  && match.VP875 <= VP875_UB  && 
//              match.HPTG1  >= HPTG1_LB  && match.HPTG1 <= HPTG1_UB  && 
//              match.VPTG1  >= VPTG1_LB  && match.VPTG1 <= VPTG1_UB  && 
//              match.HPTG2  >= HPTG2_LB  && match.HPTG2 <= HPTG2_UB  && 
//              match.VPTG2  >= VPTG2_LB  && match.VPTG2 <= VPTG2_UB  && 
//              match.BTJT2  >= BTJT2_LB  && match.BTJT2  < BTJT2_UB  &&
              match.THCURR >= THCURR_LB && match.THCURR <= THCURR_UB) {
            return true;
        }
    }
    return false;
});

const SpillCut kBNBQuality([](const caf::SRSpillProxy* sr){
    if( sr->hdr.nbnbinfo){
        GLOBAL_SPILL_INFO.clear();
        for( const auto &bnbinfo : sr->hdr.bnbinfo){ 
            fill_spill_info( sr->hdr.run, bnbinfo);  
        }
    }
    for( const auto& match: GLOBAL_SPILL_INFO) {
        if( match.event != sr->hdr.evt) continue; 
        if( match.TOR860 >= TOR860_LB && \
              match.TOR875 >= TOR875_LB && \
              match.LM875A >= LM875A_LB && \
              match.LM875B >= LM875B_LB && \
              match.LM875C >= LM875C_LB && \
              match.HP875  >= HP875_LB  && match.HP875 <= HP875_UB  && \
              match.VP875  >= VP875_LB  && match.VP875 <= VP875_UB  && \
              match.HPTG1  >= HPTG1_LB  && match.HPTG1 <= HPTG1_UB  && \
              match.VPTG1  >= VPTG1_LB  && match.VPTG1 <= VPTG1_UB  && \
              match.HPTG2  >= HPTG2_LB  && match.HPTG2 <= HPTG2_UB  && \
              match.VPTG2  >= VPTG2_LB  && match.VPTG2 <= VPTG2_UB  && \
//              match.BTJT2  >= BTJT2_LB  && match.BTJT2  < BTJT2_UB  &&
              match.THCURR >= THCURR_LB && match.THCURR <= THCURR_UB) {
            return true;
        }
    }
    return false; 
});

} // stop using namespace ana
