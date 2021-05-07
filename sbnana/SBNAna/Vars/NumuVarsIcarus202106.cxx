#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  const Var kTruthIndex = SIMPLEVAR(truth.index);
  
  const Var kPrimaryEnergy = SIMPLEVAR(truth.E);
  
  const Var kMuMaxTrack([](const caf::SRSliceProxy* slc) -> float {
      float len(-5.f);
      for (auto const& trk : slc->reco.trk)
      {
	if ( (trk.truth.p.pdg == 13 || trk.truth.p.pdg == -13) && trk.truth.p.length > len ) len = trk.truth.p.length;
      }
      return len;
    });
}
