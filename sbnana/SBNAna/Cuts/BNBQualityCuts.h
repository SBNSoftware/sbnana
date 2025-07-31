//! ////////////////////////////////////////////////////////////////////////////
//! @file: BNBQualityCuts.h                                                    
//! @author: Jacob Smith (smithja)  
//! @email: jacob.a.smith@stonybrook.edu                                           
//!                                                                           
//! @brief Header file to define BNB Quality Cuts for SBN experiments.
//! ////////////////////////////////////////////////////////////////////////////

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana {
    extern const SpillCut kTOR860Cut;
    extern const SpillCut kTOR875Cut;

    extern const SpillCut kLM875ACut;
    extern const SpillCut kLM875BCut;
    extern const SpillCut kLM875CCut;

/**
    This cut is left commented out since it is hard to implement but may be of
    future use. See the BNBQualityCuts.cxx for more information.
*/
//!    extern const SpillCut kBTJT2;

    extern const SpillCut kTHCURRCut;

    //! Figure(s) of Merit: @see getBNBFoM.cxx (and getBNBFoM2.cxx) for details
    extern const SpillCut kFoMCut; 
    extern const SpillCut kFoM2Cut;

    //! ^--- INDIVIDUAL CUTS ; COMBINATION CUTS ---v

    extern const SpillCut kBNBQualityCut_FoM1;
    extern const SpillCut kBNBQualityCut_FoM2;
}