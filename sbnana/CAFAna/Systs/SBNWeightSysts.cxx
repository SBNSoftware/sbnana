#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Systs/UniverseOracle.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  UniverseWeight::UniverseWeight(const std::string& psetName, int univIdx)
    : fPSetName(psetName), fPSetIdx(-1), fUnivIdx(univIdx)
  {
  }

  // --------------------------------------------------------------------------
  double UniverseWeight::operator()(const caf::SRSliceProxy* sr) const
  {
    if(sr->truth.index < 0) return 1;

    if(fPSetIdx == -1){
      const UniverseOracle& uo = UniverseOracle::Instance();
      fPSetIdx = uo.ParameterSetIndex(fPSetName);
    }

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = sr->truth.wgt;
    if(wgts.empty()) return 1;

    const int Nwgts = wgts[fPSetIdx].univ.size();

    static bool once = true;
    if(!once && fUnivIdx >= Nwgts){
      once = false;
      std::cout << "UniverseWeight: WARNING requesting universe " << fUnivIdx << " in parameter set " << fPSetName << " which only has size " << Nwgts << ". Will wrap-around and suppress future warnings." << std::endl;
    }

    const unsigned int unividx = fUnivIdx % Nwgts;

    return wgts[fPSetIdx].univ[unividx];
  }

  // --------------------------------------------------------------------------
  SBNWeightSyst::SBNWeightSyst(const std::string& systName)
    : ISyst(systName, systName),
      fIdx(-1)
  {
  }

  // --------------------------------------------------------------------------
  void SBNWeightSyst::Shift(double x, caf::SRSliceProxy* sr, double& weight) const
  {
    if(sr->truth.index < 0) return;

    if(fIdx == -1) fIdx = UniverseOracle::Instance().SystIndex(ShortName());

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = sr->truth.wgt;
    if(wgts.empty()) return;

    const Univs u = GetUnivs(x);

    double y = 0;
    if(u.w0 != 0) y += u.w0 * wgts[fIdx].univ[u.i0];
    if(u.w1 != 0) y += u.w1 * wgts[fIdx].univ[u.i1];

    weight *= y;
  }

  // --------------------------------------------------------------------------
  SBNWeightSyst::Univs SBNWeightSyst::GetUnivs(double x) const
  {
    auto it = fUnivs.find(x);
    if(it != fUnivs.end()) return it->second;

    Univs u;
    const UniverseOracle& uo = UniverseOracle::Instance();
    // Neighbours
    double x0, x1;
    u.i0 = uo.ClosestShiftIndex(ShortName(), x, ESide::kBelow, &x0);

    if(x0 == x){
      // Found an exact match (this is OK for small integer shifts, even with
      // floating point comparisons
      u.i1 = -1;
      u.w0 = 1;
      u.w1 = 0;
      fUnivs[x] = u;
      return u;
    }

    // Otherwise we're interpolating
    u.i1 = uo.ClosestShiftIndex(ShortName(), x, ESide::kAbove, &x1);
    // Interpolation weights
    u.w0 = (x1-x)/(x1-x0);
    u.w1 = (x-x0)/(x1-x0);

    //      std::cout << ShortName() << " " << x << " sigma, found indices " << u.i0 << " and " << u.i1 << " at " << x0 << " and " << x1 << ", will use weights " << u.w0 << " and " << u.w1 << std::endl;

    // If one of the neighbours wasn't found, we fall back to just using the
    // neighbour we did find. It would probably be better to find two
    // neighbours on the same side and extrapolate.
    if(u.i0 == -1){u.i0 = u.i1; u.w0 = u.w1 = 0.5;}
    if(u.i1 == -1){u.i1 = u.i0; u.w0 = u.w1 = 0.5;}

    fUnivs[x] = u;
    return u;
  }

  // --------------------------------------------------------------------------
  std::vector<std::string> GetSBNGenieWeightNames()
  {
    // We can't ask the UniverseOracle about this, because it doesn't get
    // properly configured until it's seen its first CAF file.
    return {"GENIEReWeight_SBN_v1_multisigma_AhtBY",
        "GENIEReWeight_SBN_v1_multisigma_BhtBY",
        "GENIEReWeight_SBN_v1_multisigma_CoulombCCQE",
        "GENIEReWeight_SBN_v1_multisigma_CV1uBY",
        "GENIEReWeight_SBN_v1_multisigma_CV2uBY",
        "GENIEReWeight_SBN_v1_multisigma_EtaNCEL",
        "GENIEReWeight_SBN_v1_multisigma_FrAbs_N",
        "GENIEReWeight_SBN_v1_multisigma_FrAbs_pi",
        "GENIEReWeight_SBN_v1_multisigma_FrCEx_N",
        "GENIEReWeight_SBN_v1_multisigma_FrCEx_pi",
        "GENIEReWeight_SBN_v1_multisigma_FrInel_N",
        "GENIEReWeight_SBN_v1_multisigma_FrInel_pi",
        "GENIEReWeight_SBN_v1_multisigma_FrPiProd_N",
        "GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi",
        "GENIEReWeight_SBN_v1_multisigma_MFP_N",
        "GENIEReWeight_SBN_v1_multisigma_MFP_pi",
        "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
        "GENIEReWeight_SBN_v1_multisigma_MaNCEL",
        "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
        "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
        "GENIEReWeight_SBN_v1_multisigma_MvNCRES",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NormCCMEC",
        "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
        "GENIEReWeight_SBN_v1_multisigma_RPA_CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE",
    };
  }

  // --------------------------------------------------------------------------
  const std::vector<const ISyst*>& GetSBNGenieWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(!ret.empty()) return ret;

    for(const std::string& name: GetSBNGenieWeightNames())
      ret.push_back(new SBNWeightSyst(name));

    return ret;
  }

  // -------------------------------------------------------------------------

  std::vector<std::string> GetSBNBoosterFluxWeightNames()
  {
    // We can't ask the UniverseOracle about this, because it doesn't get
    // properly configured until it's seen its first CAF file.
    return {"expskin_Flux",
            "horncurrent_Flux",
            "nucleoninexsec_Flux",
            "nucleonqexsec_Flux",
            "nucleontotxsec_Flux",
            "pioninexsec_Flux",
            "pionqexsec_Flux",
            "piontotxsec_Flux"};
  }

  // --------------------------------------------------------------------------
  
  const std::vector<const ISyst*>& GetSBNBoosterFluxWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(!ret.empty()) return ret;

    for(const std::string& name: GetSBNBoosterFluxWeightNames())
      ret.push_back(new SBNWeightSyst(name));

    return ret;
  }

  // --------------------------------------------------------------------------
  const std::vector<const ISyst*>& GetSBNWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(ret.empty()){
      const std::vector<const ISyst*>& g = GetSBNGenieWeightSysts();
      const std::vector<const ISyst*>& f = GetSBNBoosterFluxWeightSysts();
      ret.reserve(g.size()+f.size());
      ret.insert(ret.end(), g.begin(), g.end());
      ret.insert(ret.end(), f.begin(), f.end());
    }
    return ret;
  }
  
}
