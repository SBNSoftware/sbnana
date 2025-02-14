///////////////////////////////////////////////////////////////////////////////
// File: BNBQualityCuts.h                                                    //
// Author: Jacob Smith (smithja)                                             //
// Last edited: February 12th, 2025                                          //
//                                                                           //
// Header file to define all of the BNB Quality Cuts for the ICARUS          //
// experiment.                                                               //
///////////////////////////////////////////////////////////////////////////////

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana {
    extern const SpillCut kTOR860Cut;
    extern const SpillCut kTOR875Cut;

    extern const SpillCut kLM875ACut;
    extern const SpillCut kLM875BCut;
    extern const SpillCut kLM875CCut;

    extern const SpillCut kHP875Cut;
    extern const SpillCut kVP875Cut;

    extern const SpillCut kHPTG1Cut; 
    extern const SpillCut kVPTG1Cut;

    extern const SpillCut kHPTG2Cut;
    extern const SpillCut kVPTG2Cut;

/*
    This cut is left commented out since it is hard to implement but may be of
    future use. See the BNBQualityCuts.cxx for more information.
*/
//     extern const SpillCut kBTJT2;

    extern const SpillCut kTHCURRCut;

    // ^--- INDIVIDUAL CUTS ; COMBINATION CUTS ---v

    extern const SpillCut kBNBQualityCut_noBPMs_10Feb2025 = kTOR860Cut && kTOR875Cut && kLM875ACut && kLM875BCut && kLM875CCut && kTHCURRCut;
}
