#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana {

extern const SpillCut kTOR860;
extern const SpillCut kTOR875;

extern const SpillCut kLM875A;
extern const SpillCut kLM875B;
extern const SpillCut kLM875C;

extern const SpillCut kHP875;
extern const SpillCut kVP875;

extern const SpillCut kHPTG1; 
extern const SpillCut kVPTG1;

extern const SpillCut kHPTG2;
extern const SpillCut kVPTG2;

/*
   This cut is left commented out since it is hard to implement but may be of
   future use. See the BNBQualityCuts_<Nov2024, Jan2025>.cxx files for more
   information. -- Jacob Smith

extern const SpillCut kBTJT2; 
*/

extern const SpillCut kTHCURR;

// ^--- INDIVIDUAL CUTS ; MASTER CUT ---v

extern const SpillCut kBNBQuality_Jan25; // combination of ALL cuts given above
}
