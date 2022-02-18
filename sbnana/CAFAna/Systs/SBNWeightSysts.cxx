#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnana/CAFAna/Systs/UniverseOracle.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  UniverseWeight::UniverseWeight(const std::vector<std::string>& systs, int univIdx)
    : fNames(systs), fUnivIdx(univIdx)
  {
  }

  // --------------------------------------------------------------------------
  UniverseWeight::UniverseWeight(const std::vector<const ISyst*>& systs, int univIdx)
    : fUnivIdx(univIdx)
  {
    for(const ISyst* s: systs) fNames.push_back(s->ShortName());
  }

  // --------------------------------------------------------------------------
  double UniverseWeight::operator()(const caf::SRSliceProxy* sr) const
  {
    if(sr->truth.index < 0) return 1;

    if(fSystIdxs.empty()){
      const UniverseOracle& uo = UniverseOracle::Instance();
      for(const std::string& name: fNames){

        fSystIdxs.push_back(uo.SystIndex(name));
        fUnivOffsets.push_back(0);
        fUnivCuts.push_back(kNoCut);

      }
    }

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = sr->truth.wgt;
    if(wgts.empty()) return 1;

    // This hack improves throughput vastly
    /*
    if(fUnivIdx == 0){
      for(unsigned int i = 0; i < fNames.size(); ++i){
        for(const auto& b: wgts[fUnivIdx].univ) (void)((float)b);
      }
    }
    */

    double w = 1;

    for(unsigned int i = 0; i < fNames.size(); ++i){
      if(fUnivCuts[i](sr)){
        const unsigned int idx = fSystIdxs[i];

        // TODO: might want to "wrap around" differently in different systs to
        // avoid unwanted correlations between systs with the same number of
        // universes.
        const unsigned int unividx = (fUnivIdx + fUnivOffsets[i]) % wgts[idx].univ.size();

        w *= wgts[idx].univ[unividx];
      }
    }

    return w;
  }

  // --------------------------------------------------------------------------
  SBNWeightSyst::SBNWeightSyst(const std::string& systName,
                               const std::string& latexName,
                               const caf::ReweightType_t& type,
                               const std::string& knobName,
                               const SliceCut& cut)
    : ISyst(systName, latexName),
      fIdx(-1),
      fType(type),
      fKnobName(knobName.empty() ? systName : knobName),
      fCut(cut)
  {
  }

  // --------------------------------------------------------------------------
  void SBNWeightSyst::Shift(double x, caf::SRSliceProxy* sr, double& weight) const
  {
    if(sr->truth.index < 0) return;
    if(!fCut(sr)) return;

    if(fIdx == -1) fIdx = UniverseOracle::Instance().SystIndex(fKnobName);

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = sr->truth.wgt;
    if(wgts.empty()) return;

    const Univs u = GetUnivs(x);

    const double y0 = wgts[fIdx].univ[u.i0];
    const double y1 = wgts[fIdx].univ[u.i1];

    weight *= u.w0*y0 + u.w1*y1;
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
    u.i0 = uo.ClosestShiftIndex(fKnobName, x, ESide::kBelow, &x0);
    u.i1 = uo.ClosestShiftIndex(fKnobName, x, ESide::kAbove, &x1);
    // Interpolation weights

    //==== When the requested x exists in UniverseOracle::fShiftVals, we get x1==x2
    //==== e.g., x=1 with multisigma modes
    //==== We are comparing two doubles, so "==" is dangerous
    //==== for now using absolute difference < tolerance
    if( (fabs(x1-x0)<1e-5) ){
      u.w0 = 0.5;
      u.w1 = 0.5;
      return u;
    }

    u.w0 = (x1-x)/(x1-x0);
    u.w1 = (x-x0)/(x1-x0);

    //      std::cout << fKnobName << " " << x << " sigma, found indices " << u.i0 << " and " << u.i1 << " at " << x0 << " and " << x1 << ", will use weights " << u.w0 << " and " << u.w1 << std::endl;

    // If one of the neighbours wasn't found, we fall back to just using the
    // neighbour we did find. It would probably be better to find two
    // neighbours on the same side and extrapolate.
    if(u.i0 == -1){u.i0 = u.i1; u.w0 = u.w1 = 0.5;}
    if(u.i1 == -1){u.i1 = u.i0; u.w0 = u.w1 = 0.5;}

    fUnivs[x] = u;
    return u;
  }

  // --------------------------------------------------------------------------

  std::vector<std::string> GetSBNGenieWeightPSet(const caf::ReweightType_t& rwType){

    std::vector<std::string> names_multisim = {
    "sbnd",
    };

    std::vector<std::string> names_multisigma = {
    "AhtBY",
    "BhtBY",
    "CV1uBY",
    "CV2uBY",
    "EtaNCEL",
    "FormZone",
    "FrAbs_N",
    "FrAbs_pi",
    "FrCEx_N",
    "FrCEx_pi",
    "FrInel_N",
    "FrInel_pi",
    "FrPiProd_N",
    "FrPiProd_pi",
    "MFP_N",
    "MFP_pi",
    "MaCCQE",
    "MaCCRES",
    "MaNCEL",
    "MaNCRES",
    "MvCCRES",
    "MvNCRES",
    "NonRESBGvbarnCC1pi",
    "NonRESBGvbarnCC2pi",
    "NonRESBGvbarnNC1pi",
    "NonRESBGvbarnNC2pi",
    "NonRESBGvbarpCC1pi",
    "NonRESBGvbarpCC2pi",
    "NonRESBGvbarpNC1pi",
    "NonRESBGvbarpNC2pi",
    "NonRESBGvnCC1pi",
    "NonRESBGvnCC2pi",
    "NonRESBGvnNC1pi",
    "NonRESBGvnNC2pi",
    "NonRESBGvpCC1pi",
    "NonRESBGvpCC2pi",
    "NonRESBGvpNC1pi",
    "NonRESBGvpNC2pi",
    };

    if(rwType==caf::kMultiSim) return names_multisim;
    else if(rwType==caf::kMultiSigma) return names_multisigma;
    else{
      return std::vector<std::string>();
    }


  }

  std::vector<const ISyst*> GetSBNGenieWeightSysts(const caf::ReweightType_t& rwType)
  {
    std::vector<const ISyst*> ret;
    if(!ret.empty()) return ret;

    // We can't ask the UniverseOracle about this, because it doesn't get
    // properly configured until it's seen its first CAF file.

    const std::vector<std::string> names = GetSBNGenieWeightPSet(rwType);

    for(const std::string& name: names){
      std::string psetname(name);
      if(rwType==caf::kMultiSim){
        psetname = "genie_"+name+"_multisim_Genie";
      }
      else if(rwType==caf::kMultiSigma){
        psetname = "genie_"+name+"_multisigma_Genie";
      }
      ret.push_back(new SBNWeightSyst(psetname, name, rwType));
    }

    return ret;
  }

  // --------------------------------------------------------------------------
  /*
  const std::vector<const ISyst*>& GetBoosterFluxWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(ret.empty()){
      const UniverseOracle& uo = UniverseOracle::Instance();
      ret.reserve(uo.Systs().size());
      for(const std::string& name: uo.Systs()){
        if(name.find("genie") != std::string::npos) continue;
        ret.push_back(new SBNWeightSyst(name));
      }
    }
    return ret;
  }

  // --------------------------------------------------------------------------
  const std::vector<const ISyst*>& GetSBNWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(ret.empty()){
      const std::vector<const ISyst*>& g = GetSBNGenieWeightSysts();
      const std::vector<const ISyst*>& f = GetBoosterFluxWeightSysts();
      ret.reserve(g.size()+f.size());
      ret.insert(ret.end(), g.begin(), g.end());
      ret.insert(ret.end(), f.begin(), f.end());
    }
    return ret;
  }
  */
}
