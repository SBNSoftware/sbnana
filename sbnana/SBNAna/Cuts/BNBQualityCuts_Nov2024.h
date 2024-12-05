#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana {

extern const SpillCut kPrelimTOR860; 
extern const SpillCut kPrelimTOR875;

extern const SpillCut kPrelimLM875A;
extern const SpillCut kPrelimLM875B;
extern const SpillCut kPrelimLM875C;

extern const SpillCut kPrelimTHCURR;


extern const SpillCut kPRELIMINARY; // combination of all preliminary cuts

// ^--- PRELIMINARY CUTS ; PHYSICALLY MEANINGFUL CUTS ---v

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

extern const SpillCut kBTJT2; 

extern const SpillCut kTHCURR;

extern const SpillCut kPHYSICALLY_MEANINGFUL; // combination of all physically 
                                     // meaningful cuts

// ^--- PHYSICALLY MEANINGFUL CUTS ; MASTER CUT ---v

extern const SpillCut kBNBQuality_Nov24; // combination of ALL cuts given above
}
