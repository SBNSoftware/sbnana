#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Vars/Vars.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <cmath>

namespace ana
{
    const Var kHasTruthMatch(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( slc->truth.index != -1);
      }
      );

    const Var kCompletness(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)slc->tmatch.eff : -5.f);
      }
      );

    const Var kTruthEnergy(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.E : -5.f );
      }
      );

    const Var kTruthVtxX(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.position.x : -999.f );
      }
      );

    const Var kTruthVtxY(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.position.y : -999.f );
      }
      );

    const Var kTruthVtxZ(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)slc->truth.position.z : -999.f );
      }
      );

    const Var kTruthVtxDistX(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)std::abs( slc->truth.position.x - kSlcVtxX(slc) ) : -999.f );
      }
      );

    const Var kTruthVtxDistY(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)std::abs( slc->truth.position.y - kSlcVtxY(slc) ) : -999.f );
      }
      );

    const Var kTruthVtxDistZ(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)std::abs( slc->truth.position.z - kSlcVtxZ(slc) ) : -999.f );
      }
      );

    const Var kTruthVtxDistMag(
      [](const caf::SRSliceProxy *slc) -> double
      {
        return ( kHasTruthMatch(slc) ? (float)std::hypot( kTruthVtxDistX(slc), kTruthVtxDistY(slc), kTruthVtxDistZ(slc) ) : -5.f );
      }
      );

    // -----------------------------------------------------------------
    // Spill vars

    const SpillVar kTruthNuEnergy([](const caf::SRSpillProxy* sr)  -> double {
      return (sr->mc.nnu != 1 ? -5.f : (float)sr->mc.nu[0].E);
    });

    const SpillVar kTruthLeptonEnergy([](const caf::SRSpillProxy* sr)  -> double {
      if (sr->mc.nnu != 1) return -5.f;

      for (auto const& prim : sr->mc.nu[0].prim) {
        if (std::abs(prim.pdg) != 11 && std::abs(prim.pdg) != 13)
          continue;

        return prim.startE;
      }
      return -5.f;
    });

}
