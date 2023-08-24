#include "sbnana/SBNAna/Cuts/NuMIXSecTrueEventVars.h"

using namespace ana::PrimaryUtil;
using namespace caf;

namespace ana{

  std::vector<double> GetTrueVarVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const SRTrueInteractionProxy&)> isSignal,
    std::function<double(const SRTrueInteractionProxy&) > trueth_var
  ){

    std::vector<double> vals;

    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi
      vals.push_back( trueth_var(nu) );
    }

    return vals;

  }

  std::vector<double> GetCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const SRTrueInteractionProxy&)> isSignal
  ){

    vector<double> vals;

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
    std::function<bool(const SRTrueInteractionProxy&)> isSignal
  ){

    vector<double> vals;
    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi
      vals.push_back( FluxWeightNuMI.GetNuWeight(nu) );
    }
    return vals;

  }

}
