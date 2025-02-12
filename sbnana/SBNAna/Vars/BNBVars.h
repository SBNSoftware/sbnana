///////////////////////////////////////////////////////////////////////////////
// File: BNBVars.h                                                           //
// Author: Jacob Smith (smithja)                                             //
// Last edited: February 12th, 2025                                          //
//                                                                           //
// Header file to define all of the variables used in BNB Quality Cuts for   //
// the ICARUS experiment.                                                    //
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{
    // readout time for a given spill
    // units: s (i.e. seconds)
    extern const SpillVar kSpillTimeSec;

    // toroid devices: monitor protons on target (POT)
    // units: POT
    extern const SpillVar kSpillTOR860;
    extern const SpillVar kSpillTOR875;
    
    // loss monitors: measure backscattering of beam just upstream of target
    // NB: larger reading is better, implies beam is hitting more of target
    // units: rad / s (i.e. rads per second)
    extern const SpillVar kSpillLM875A;
    extern const SpillVar kSpillLM875B;
    extern const SpillVar kSpillLM875C;
    
    // beam position monitors: measure beam position at different points upstream
    // of target
    // NB: not always centered at (0, 0); check SBNAna/Cuts/BNBQualityCuts.cxx
    //     for information on nominal beam position for different runs
    // units: mm
    extern const SpillVar kSpillHP875;
    extern const SpillVar kSpillVP875;
    
    extern const SpillVar kSpillHPTG1;
    extern const SpillVar kSpillVPTG1;
    
    extern const SpillVar kSpillHPTG2;
    extern const SpillVar kSpillVPTG2;
    
    // DEPRECIATED, SEE NB
    // temperature reading near target
    // NB: this beam monitoring device provided diminishing returns when combined
    //     with all of the other devices listed here. additionally, temperature
    //     readings can be correlated between spills if they are close enough in
    //     time
    // units: degrees celsius
//    extern const SpillVar kSpillBTJT2;
    
    // focusing horn current
    // units: kA
    extern const SpillVar kSpillTHCURR;
}
