#include "sbnana/SBNAna/Cuts/TrueParticleHelpers.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana{

  bool IsPrimary(const caf::SRTrueParticleProxy& p) {
    return p.start_process == caf::kG4primary;
  }

  bool HasBraggPeak(const caf::SRTrueParticleProxy& p) {
    // check contained, id, & end process (stopping-only)
    return p.contained &&
      ( abs(p.pdg) == 13 || abs(p.pdg) == 2212 ||
        abs(p.pdg) == 211 || abs(p.pdg) == 321 ) &&
      ( p.end_process ==  caf::kG4CoupledTransportation ||
        p.end_process ==  caf::kG4FastScintillation ||
        p.end_process ==  caf::kG4Decay ||
        p.end_process ==  caf::kG4muMinusCaptureAtRest );
  }

  bool IsGenie(const caf::SRTrueParticleProxy& p) {
    return p.gstatus !=  caf::kNotGenie;
  }

  bool IsStable(const caf::SRTrueParticleProxy& p) {
    // non-genie particles are stable, otherwise check with genie
    return ( p.gstatus ==  caf::kNotGenie || p.gstatus ==  caf::kIStStableFinalState);
  }

}
