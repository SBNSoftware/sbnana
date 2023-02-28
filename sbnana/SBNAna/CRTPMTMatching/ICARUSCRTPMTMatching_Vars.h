#include "sbnana/SBNAna/CRTPMTMatching/ICARUSCRTPMTMatching_Tool.h"

using namespace std;
using namespace ana;

namespace ICARUSCRTPMTMatching{

  // - PMT-CRT matching
  extern const SpillMultiVar spillvarOpFlashTime;
  extern const SpillMultiVar spillvarValidOpFlashTime;
  extern const SpillMultiVar spillvarInTimeOpFlashTime; // Valid AND InTime
  extern const SpillMultiVar spillvarCRTPMTTime;
  extern const SpillMultiVar spillvarCRTPMTMatchingID;

  extern const SpillMultiVar spillvarMatchID2_CRTHitPosXs;

} // end namespace ICARUSCRTPMTMatching
