#pragma once

#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"

#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSCRTPMTMatching{

  // TODO Use spill info in the future..
  enum GateType
  {
    UNKNOWN = 0,
    BNB = +1,
    OffBeamBNB = -1,
    NUMI = +2,
    OffBeamNUMI = -2,
  };

  class CRTPMTMatchingTool{

    public:

    CRTPMTMatchingTool();

    static CRTPMTMatchingTool& Instance();

    void SetGateType(GateType gt) const;
    void SetInTimeRange(double t_min, double t_max) const;
    bool IsInTime(double t_gate) const;
    std::vector<int> GetMatchedCRTHitIndex(
      double opt,
      const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
      int mode=0
    ) const;
    int GetMatchID(
      double opt,
      const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
      int mode=0
    ) const;

    bool IsNegativeTOF(double timediff) const;

    mutable bool Debug;

    mutable bool UseTS0;
    mutable GateType GT;
    mutable double timecut_min, timecut_max;

  };

  // instance
  static const CRTPMTMatchingTool& cpmt = CRTPMTMatchingTool::Instance();

}


