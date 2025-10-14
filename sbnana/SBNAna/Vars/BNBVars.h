//! ////////////////////////////////////////////////////////////////////////////
//! @file: BNBVars.h                                                           
//! @author: Jacob Smith (smithja)
//! @email: jacob.a.smith@stonybrook.edu                                       
//!                                                                           
//! @brief Header file to define variables used in BNB Quality Cuts for the 
//! SBN experiments.                                                    
//! ////////////////////////////////////////////////////////////////////////////

#pragma once

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{
    //! Readout time for a given spill.
    //! units: s (i.e. seconds)
    extern const SpillVar kSpillTimeSec;

    //! Toroid Devices: monitor protons on target (POT)
    //! units: POT
    extern const SpillVar kSpillTOR860;
    extern const SpillVar kSpillTOR875;
    extern const SpillVar kSpillTOR;
    
    //! Loss Monitors: measure backscattering of beam just upstream of target
    //! @note larger reading is better, implies beam is hitting more of target
    //! units: rad / s (i.e. rads per second)
    extern const SpillVar kSpillLM875A;
    extern const SpillVar kSpillLM875B;
    extern const SpillVar kSpillLM875C;
    extern const SpillVar kSpillLM875;

    //! Beam Position Monitors: measure beam position at different points
    //! upstream of target
    //! @note not always centered at (0, 0); @see SBNAna/Cuts/BNBQualityCuts.cxx
    //! for information on nominal beam position for different runs
    //! @note VPTG1 and VPTG2 unreliable for ICARUS (and presumably SBND) runs;
    //! VP873 identified as a suitable replacement
    //! units: mm
    extern const SpillVar kSpillHP875;
    extern const SpillVar kSpillVP875;
    
    extern const SpillVar kSpillHPTG1;
 //!   extern const SpillVar kSpillVPTG1;
    
    extern const SpillVar kSpillHPTG2;
 //!   extern const SpillVar kSpillVPTG2;

    extern const SpillVar kSpillVP873;

    //! @deprecated
    //! Temperature reading near target.
    //! @note this beam monitoring device provided diminishing returns when 
    //! combined with all of the other devices listed here. additionally, 
    //! temperature readings can be correlated between spills if they are close 
    //! enough in time
    //! units: degrees celsius
//!    extern const SpillVar kSpillBTJT2;
    
    //! Focusing Horn Current
    //! units: kA
    extern const SpillVar kSpillTHCURR;

    //! Multi-wire Readout Fit Parameters:
    //! @note We determine beam sigma primarily by fitting gaussian curves to
    //! multi-wire device data; we can also extract beam position from fits.
    //! Beam position should be taken from beam position monitor variables
    //! above, but can be taken from below if needed.
    //! @note Sigmas and positions provided for both horizontal and vertical
    //! directions perpendicular to the beam direction.
    //! units: mm
    extern const SpillVar kSpillMW875HorSigma;
    extern const SpillVar kSpillMW875VerSigma;
    extern const SpillVar kSpillMW876HorSigma;
    extern const SpillVar kSpillMW876VerSigma;

    //! Figure(s) of Merit: @see getBNBFoM.cxx for details
    extern const SpillVar kSpillFoM_noMultiWire;
    extern const SpillVar kSpillFoM;
}