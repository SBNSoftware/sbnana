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

  std::vector<double> GetSpillCutTypeVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal
  ){

    std::vector<double> vals;

    bool ThisSpillPass = kNuMIValidTrigger(sr);

    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi

      if(ThisSpillPass) vals.push_back( 1 );
      else vals.push_back( 0 );

    }

    return vals;

  }

  std::vector<double> GetCutTypeWithoutShowerCutVectorPerNu(
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

        int this_slice_cuttpye = kNuMICutTypeWithoutShowerCut(&slc);
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

  std::vector<double> GetSigmaWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    const std::string& psetName,
    double shift // in sigma
  ){

    std::vector<double> vals;
    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi

      const UniverseOracle& uo = UniverseOracle::Instance();
      int fPSetIdx = uo.ParameterSetIndex(psetName);

      const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = nu.wgt;

      double this_weight(1.);
      if(!wgts.empty()){
        // Neighbours
        double x0, x1;
        Univs u;
        u.i0 = uo.ClosestShiftIndex(psetName, shift, ESide::kBelow, &x0);
        if(x0 == shift){
          // Found an exact match (this is OK for small integer shifts, even with
          // floating point comparisons
          u.i1 = -1;
          u.w0 = 1;
          u.w1 = 0;
        }
        else{
          // Otherwise we're interpolating
          u.i1 = uo.ClosestShiftIndex(psetName, shift, ESide::kAbove, &x1);
          // Interpolation weights
          u.w0 = (x1-shift)/(x1-x0);
          u.w1 = (shift-x0)/(x1-x0);
        }

        // If one of the neighbours wasn't found, we fall back to just using the
        // neighbour we did find. It would probably be better to find two
        // neighbours on the same side and extrapolate.
        if(u.i0 == -1){u.i0 = u.i1; u.w0 = u.w1 = 0.5;}
        if(u.i1 == -1){u.i1 = u.i0; u.w0 = u.w1 = 0.5;}

        double y = 0; 
        if(u.w0 != 0) y += u.w0 * wgts[fPSetIdx].univ[u.i0];
        if(u.w1 != 0) y += u.w1 * wgts[fPSetIdx].univ[u.i1];
        this_weight = y;
      }
      vals.push_back( this_weight );
    }
    return vals;

  }
  SpillMultiVar GetSigmaWeightSpillMultiVarPerSignalNu(
    const std::string& psetName,
    double shift // in sigma
  ){
    return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetSigmaWeightVectorPerNu(sr, Is1muNp0piWithPhaseSpaceCut, psetName, shift);
    });
  }

  std::vector<double> GetUniverseWeightVectorPerNu(
    const caf::SRSpillProxy* sr,
    std::function<bool(const caf::SRTrueInteractionProxy&)> isSignal,
    const std::string& psetName,
    int univIdx
  ){

    std::vector<double> vals;
    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi

      const UniverseOracle& uo = UniverseOracle::Instance();
      int fPSetIdx = uo.ParameterSetIndex(psetName);

      const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = nu.wgt;

      double this_weight(1.);
      if(!wgts.empty()){
        const int Nwgts = wgts[fPSetIdx].univ.size();
        const unsigned int unividx = univIdx % Nwgts;
        this_weight = wgts[fPSetIdx].univ[unividx];
      }
      vals.push_back( this_weight );
    }
    return vals;

  }
  SpillMultiVar GetUniverseWeightSpillMultiVarPerSignalNu(
    const std::string& psetName,
    int univIdx
  ){
    return SpillMultiVar( [=](const caf::SRSpillProxy *sr) -> std::vector<double> {
      return GetUniverseWeightVectorPerNu(sr, Is1muNp0piWithPhaseSpaceCut, psetName, univIdx);
    });
  }

}
