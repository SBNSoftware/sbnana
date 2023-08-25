#include "sbnana/SBNAna/Cuts/NuMIXSecTrueEventVars.h"

namespace ana{

  std::vector<double> GetTrueVarVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    std::function<double(const caf::SRTrueInteractionProxy&) > TrueVar
  ){

    std::vector<double> vals;

    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi
      vals.push_back( TrueVar(nu) );
    }

    return vals;

  }

  SpillMultiVar GetTrueSpillMultiVarPerSignalNu(
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    std::function<double(const caf::SRTrueInteractionProxy&) > TrueVar
  ){

    return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetTrueVarVectorPerNu(sr, isSignal, TrueVar);
    });

  }

  std::vector<double> GetCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  ){

    std::vector<double> vals;

    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi

      int this_nu_index = nu.index;

      bool ThisNuPassSignal = false;

      // Loop over reco slice
      for ( auto const& slc : sr->slc ) {
        if ( slc.truth.index < 0 ) continue;
        else if ( slc.truth.index != this_nu_index ) continue;

        int this_slice_cuttpye = kNuMICutType(&slc);
        if(this_slice_cuttpye==1){
          ThisNuPassSignal = true;
          break;
        }
      }

      if(ThisNuPassSignal) vals.push_back( 1 );
      else vals.push_back( 0 );

    }

    return vals;

  }

  std::vector<double> GetNuMIPPFXWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  ){

    std::vector<double> vals;
    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi
      vals.push_back( FluxWeightNuMI.GetNuWeight(nu) );
    }
    return vals;

  }

}
