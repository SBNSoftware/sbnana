///////////////////////////////////////////////////////////////////////////////
// File: BNBVars.cxx                                                         //
// Author: Jacob Smith (smithja)                                             //
// Last edited: February 12th, 2025                                          //
//                                                                           //
// Explicit definitions of the variables used in BNB Quality Cuts for the    //
// ICARUS experiment. These SpillVars are declared with Bruce Howard's       //
// SIMPLESPILLVAR() function allowing SpillCuts to be easily integrated with //
// creating Spectra (like what is done with regular Cuts).                   //
///////////////////////////////////////////////////////////////////////////////

#include "sbnana/SBNAna/Vars/BNBVars.h"

#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <iostream>

namespace ana
{
    // readout time for a given spill
    // units: s (i.e. seconds)
    const SpillVar kSpillTimeSec = SIMPLESPILLVAR( hdr.spillbnbinfo.spill_time_sec);

    // toroid devices: monitor protons on target (POT)
    // units: POT
    const SpillVar kSpillTOR860 = SIMPLESPILLVAR( hdr.spillbnbinfo.TOR860);
    const SpillVar kSpillTOR875 = SIMPLESPILLVAR( hdr.spillbnbinfo.TOR875);
    
    // loss monitors: measure backscattering of beam just upstream of target
    // NB: larger reading is better, implies beam is hitting more of target
    // units: rad / s (i.e. rads per second)
    const SpillVar kSpillLM875A = SIMPLESPILLVAR( hdr.spillbnbinfo.LM875A);
    const SpillVar kSpillLM875B = SIMPLESPILLVAR( hdr.spillbnbinfo.LM875B);
    const SpillVar kSpillLM875C = SIMPLESPILLVAR( hdr.spillbnbinfo.LM875C);
    
    // beam position monitors: measure beam position at different points upstream
    // of target
    // NB: not always centered at (0, 0); check SBNAna/Cuts/BNBQualityCuts.cxx
    //     for information on nominal beam position for different runs
    // units: mm
    const SpillVar kSpillHP875 = SIMPLESPILLVAR( hdr.spillbnbinfo.HP875);
    const SpillVar kSpillVP875 = SIMPLESPILLVAR( hdr.spillbnbinfo.VP875);
   
    const SpillVar kSpillHPTG1 = SIMPLESPILLVAR( hdr.spillbnbinfo.HPTG1);
    const SpillVar kSpillVPTG1 = SIMPLESPILLVAR( hdr.spillbnbinfo.VPTG1);

    const SpillVar kSpillHPTG2 = SIMPLESPILLVAR( hdr.spillbnbinfo.HPTG2);
    const SpillVar kSpillVPTG2 = SIMPLESPILLVAR( hdr.spillbnbinfo.VPTG2);
 
    // DEPRECIATED, SEE NB
    // temperature reading near target
    // NB: this beam monitoring device provided diminishing returns when combined
    //     with all of the other devices listed here. additionally, temperature
    //     readings can be correlated between spills if they are close enough in
    //     time
    // units: degrees celsius
//    const SpillVar kSpillBTJT2 = SIMPLESPILLVAR( hdr.spillbnbinfo.BTJT2);
    
    // focusing horn current
    // units: kA
    const SpillVar kSpillTHCURR = SIMPLESPILLVAR( hdr.spillbnbinfo.THCURR);
}
