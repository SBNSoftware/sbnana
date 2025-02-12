///////////////////////////////////////////////////////////////////////////////
// File: BNBQualityCuts.h                                                    //
// Author: Jacob Smith (smithja)                                             //
// Last edited: March 6th, 2025                                              //
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

/*
    This cut is left commented out since it is hard to implement but may be of
    future use. See the BNBQualityCuts.cxx for more information.
*/
//     extern const SpillCut kBTJT2;

    extern const SpillCut kTHCURRCut;

    // ^--- INDIVIDUAL CUTS ; COMBINATION CUTS ---v

    extern const SpillCut kBNBQualityCut;
}
