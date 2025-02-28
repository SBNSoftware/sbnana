///////////////////////////////////////////////////////////////////////////////
// File: BNBQualityCuts.cxx                                                  //
// Author: Jacob Smith (smithja)                                             //
// Last edited: February 12th, 2025                                          //
//                                                                           //
// Explicit definitions of the BNB Quality Cuts used for the ICARUS          //
// experiment.                                                               //
///////////////////////////////////////////////////////////////////////////////

#include "sbnana/CAFAna/Core/Cut.h"

#include "sbnana/SBNAna/Vars/BNBVars.h"
#include "sbnana/SBNAna/Cuts/BNBQualityCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana { // use namespace ana for entire file
/*
    While there is sometimes discernable structure in a given BNB quality 
    variable (e.g. TOR860 and TOR875 stay around 4.5E12 POT when only the BNB is
    on and 3.0E12 POT when both the BNB and NuMI beams are on), exact behavior 
    of such variables is not known in detail as a function of time. Thus, the 
    latest (read 'final') version of a given cut will be the most stringent cut
     we can make with a limited knowledge of the detailed behavior of that cut.
*/
const double TOR860_LB = +100e9; // units: POT
const double TOR875_LB = +100e9; 

const double LM875A_LB = +1e-2; // units: rad/s
const double LM875B_LB = +1e-2;
const double LM875C_LB = +1e-2;

/* 
    SMITHJA: I am waiting on the External Beams group to provide a table
    of the nominal beam position monitor (BPM) variables for all of the ICARUS
    data runs. Until that time, the cuts for the BPM variables are the nominal
    BPM values from MicroBooNE's getFOM2.cxx file created in part by Zarko
    Pavlovic and located at 
    https://github.com/uboone/ubraw/blob/develop/ubraw/BeamData/getFOM2.cxx 

    Upper and lower bound windows can be tweaked as needed. But if we continue
    to use the nominal BPM values from MicroBooNE, these should be encoded in
    the offsets below. These offsets are taken from getFOM2.cxx linked above.
*/
double hp875_offset = -3.40271; // units: mm
double vp875_offset = +1.47511;
double hptg1_offset = +0.457014;
double vptg1_offset = +0.38914;
double hptg2_offset = +0.8;
double vptg2_offset = +1.0;

/* 
    SMITHJA: The MicroBooNE offsets allow no spills to pass through. We are
    instead re-setting the offsets to zero and hoping External Beams will 
    come through...
*/
hp875_offset = 0; // units: mm
vp875_offset = 0;
hptg1_offset = 0;
vptg1_offset = 0;
hptg2_offset = 0;
vptg2_offset = 0;

const double HP875_LB = -0.25, HP875_UB = +0.25; // units: mm
const double VP875_LB = -0.25, VP875_UB = +0.30;
const double HPTG1_LB = -0.25, HPTG1_UB = +0.25;
const double VPTG1_LB = -0.20, VPTG1_UB = +0.30;
const double HPTG2_LB = -0.25, HPTG2_UB = +0.25;
const double VPTG2_LB = -0.20, VPTG2_UB = +0.30;

/*
    The BTJT2 variable is somewhat unreliable to use as a spill-by-spill cut.
    BTJT2 measures the temperature of the BNB near the target. Since runs can be
    started/stopped more quickly than the thermal fluctuations of BTJT2, there
    can be some correlation between spills. For this reason BTJT2 is deemed to 
    have diminishing returns compared to other BNB quality variables listed. 
    You may employ cuts on BTJT2 at your own peril.
*/

// const double BTJT2_LB = 999, BTJT2_UB = 999; // units: deg C


/* 
    SMITHJA: Further cuts on THCURR could/should be made since it was discovered 
    that the horn current is not always centered around 174 kA. 
    Further work is needed to understand horn current for different runs. The
    cut ranging from 173 to 175 kA is for ICARUS Run 2 data. Generally it was
    found that a range of 173 to 177 included the horn currents for a majority
    of ICARUS data runs.
*/
// const double THCURR_LB = +173, THCURR_UB = +177; // units: kA
const double THCURR_LB = +173, THCURR_UB = +175; // units: kA

//-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

const SpillCut kTOR860Cut = kSpillTOR860 >= TOR860_LB;
const SpillCut kTOR875Cut = kSpillTOR875 >= TOR875_LB;

const SpillCut kLM875ACut = kSpillLM875A >= LM875A_LB;
const SpillCut kLM875BCut = kSpillLM875B >= LM875B_LB;
const SpillCut kLM875CCut = kSpillLM875C >= LM875C_LB;

const SpillCut kHP875Cut  ( [](const caf::SRSpillProxy* sr) { 
    return ( kSpillHP875(sr) + hp875_offset) >= HP875_LB && ( kSpillHP875(sr) + hp875_offset) <= HP875_UB;
} );

const SpillCut kVP875Cut  ( [](const caf::SRSpillProxy* sr) { 
    return ( kSpillVP875(sr) + vp875_offset) >= VP875_LB && ( kSpillVP875(sr) + vp875_offset) <= VP875_UB;
} );

const SpillCut kHPTG1Cut  ( [](const caf::SRSpillProxy* sr) { 
    return ( kSpillHPTG1(sr) + hptg1_offset) >= HPTG1_LB && ( kSpillHPTG1(sr) + hptg1_offset) <= HPTG1_UB;
} );

const SpillCut kVPTG1Cut  ( [](const caf::SRSpillProxy* sr) { 
    return ( kSpillVPTG1(sr) + vptg1_offset) >= VPTG1_LB && ( kSpillVPTG1(sr) + vptg1_offset) <= VPTG1_UB;
} );

const SpillCut kHPTG2Cut  ( [](const caf::SRSpillProxy* sr) { 
    return ( kSpillHPTG2(sr) + hptg2_offset) >= HPTG2_LB && ( kSpillHPTG2(sr) + hptg2_offset) <= HPTG2_UB;
} );

const SpillCut kVPTG2Cut  ( [](const caf::SRSpillProxy* sr) { 
    return ( kSpillVPTG2(sr) + vptg2_offset) >= VPTG2_LB && ( kSpillVPTG2(sr) + vptg2_offset) <= VPTG2_UB;
} );

// DEPRECIATED: see above notes about BTJT2
//const SpillCut kBTJT2Cut = kSpillBTJT2 >= BTJT2_LB && kSpillBTJT2 <= BTJT2_UB;

const SpillCut kTHCURRCut = kSpillTHCURR >= THCURR_LB && kSpillTHCURR <= THCURR_UB;

// ^--- INDIVIDUAL CUTS ; COMBINATION CUTS ---v

const SpillCut kBNBQualityCut_noBPMs_10Feb2025 = kTOR860Cut && kTOR875Cut && kLM875ACut && kLM875BCut && kLM875CCut && kTHCURRCut;

// WARNING: Beam position monitor (BPM) cuts are not finalized as of February
//          20th, 2025. 
const SpillCut kBNBQualityCut_20Feb2025 = kBNBQualityCut_noBPMs_10Feb2025 && kHP875Cut && kVP875Cut && kHPTG1Cut && kVPTG1Cut && kHPTG2Cut && kVPTG2Cut;

} // stop using namespace ana
