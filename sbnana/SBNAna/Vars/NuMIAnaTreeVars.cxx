#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Vars/NuMIAnaTreeVars.h"

#include "TVector3.h"

namespace ana {

  const Var kNuMICutType([](const caf::SRSliceProxy* slc) -> int {
    // Check if this is signal or sideband or neither
    std::vector<double> chargedpion_indices = kNuMIChargedPionCandidateIdxs(slc);
    std::vector<double> neutralpion_indices = kNuMIPhotonCandidateIdxs(slc);

    if ( chargedpion_indices.size()==0 && neutralpion_indices.size()==0 ) return 1;
    else if ( chargedpion_indices.size()>0 && neutralpion_indices.size()==0 ) return 2;
    else if ( chargedpion_indices.size()==0 && neutralpion_indices.size()>0 ) return 3;
    return 0;
  });

  const Var kNuMISliceCategory([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return 4; // Is cosmic/unmatched

    if ( !slc->truth.iscc ) return 3;

    bool failsBasic = false;
    if ( abs(slc->truth.pdg) != 14 ||
				 !slc->truth.iscc ||
				 std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z) ||
				 !isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) )
      failsBasic = true; // not signal

    // primaries:
    unsigned int nMu(0), nP(0), nPi(0);
    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      if ( abs(prim.pdg) == 13 && momentum > 0.226 ) nMu+=1;
      if ( abs(prim.pdg) == 2212 && isContainedVol(prim.end.x,prim.end.y,prim.end.z) && (momentum > 0.4 && momentum < 1.) ) nP+=1;
      if ( abs(prim.pdg) == 111 || abs(prim.pdg) == 211 ) nPi+=1;
    }

    if ( failsBasic==false && nMu==1 && nP>0 && nPi==0 ) return 1;
    else return 2;
  });

  const Var kNuMITruePDG([](const caf::SRSliceProxy* slc) -> int {
      if ( slc->truth.index < 0 ) return 0;
      return slc->truth.pdg;
    });

  // Can make this read from the CAF but...
  const Var kNuMIIsFHC([](const caf::SRSliceProxy* slc) -> int {
      return 1;
    });

  // See https://icarus-exp.fnal.gov/at_work/software/doc/icaruscode/latest/SREnums_8h_source.html#l00062
  const Var kNuMITrueMode([](const caf::SRSliceProxy* slc) -> int {
      if ( slc->truth.index < 0 ) return -1;
      return slc->truth.genie_mode;
    });

  const Var kNuMITrueTarget([](const caf::SRSliceProxy* slc) -> int {
      if ( slc->truth.index < 0 ) return 0;
      return slc->truth.targetPDG;
    });

  const Var kNuMITrueIsCC([](const caf::SRSliceProxy* slc) -> int {
      if ( slc->truth.index < 0 ) return 0;
      return slc->truth.iscc;
    });

  const Var kNuMITrueENu([](const caf::SRSliceProxy* slc) -> float {
      if ( slc->truth.index < 0 ) return -5.f;
      return slc->truth.E;
    });

}
