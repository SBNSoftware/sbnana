#include "sbnana/CAFAna/Systs/SBNOnOffSysts.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Systs/UniverseOracle.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  SBNOnOffSyst::SBNOnOffSyst(const std::string& systName)
    : ISyst(systName, systName),
      fIdx(-1)
  {
  }

  // --------------------------------------------------------------------------
  void SBNOnOffSyst::Shift(double x, caf::SRSliceProxy* sr, double& weight) const
  {
    if(sr->truth.index < 0) return;

    if(fIdx == -1) fIdx = UniverseOracle::Instance().SystIndex(ShortName());

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = sr->truth.wgt;
    if(wgts.empty()) return;

    const UniverseOracle& uo = UniverseOracle::Instance();
    double x1;
    int i = uo.ClosestShiftIndex(ShortName(), 1, ESide::kBelow, &x1);
    double wgt = wgts[fIdx].univ[i];

    if(x < 0) x *= -1;
    if(x != 0) weight *= x < 1.0 ? x * wgt : wgt;
  }

  // --------------------------------------------------------------------------
  std::vector<std::string> GetSBNOnOffNames()
  {
    // We can't ask the UniverseOracle about this, because it doesn't get
    // properly configured until it's seen its first CAF file.
    return { "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
             "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
             "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
             "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad"
    };
  }

  // --------------------------------------------------------------------------
  const std::vector<const ISyst*>& GetSBNOnOffSysts()
  {
    static std::vector<const ISyst*> ret;
    if(!ret.empty()) return ret;

    for(const std::string& name: GetSBNOnOffNames())
      ret.push_back(new SBNOnOffSyst(name));

    return ret;
  }
}
