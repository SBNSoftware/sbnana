#include <iostream>
#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana{

  const Cut kHasPrimaryMuonTrk = kPrimaryMuonTrkIdx != -1;

}
