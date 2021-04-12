#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

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
      for(const std::string& name: fNames) fSystIdxs.push_back(uo.SystIndex(name));
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
      const unsigned int idx = fSystIdxs[i];

      // TODO: might want to "wrap around" differently in different systs to
      // avoid unwanted correlations between systs with the same number of
      // universes.
      const unsigned int unividx = fUnivIdx % wgts[idx].univ.size();

      w *= wgts[idx].univ[unividx];
    }

    return w;
  }

  // --------------------------------------------------------------------------
  SBNWeightSyst::SBNWeightSyst(const std::string& name, const SliceCut& cut)
    : ISyst(name, name), fIdx(-1), fCut(cut)
  {
    //    assert(UniverseOracle::Instance().SystExists(name));
  }

  // --------------------------------------------------------------------------
  void SBNWeightSyst::Shift(double x, caf::SRSliceProxy* sr, double& weight) const
  {
    if(sr->truth.index < 0) return;
    if(!fCut(sr)) return;

    if(fIdx == -1) fIdx = UniverseOracle::Instance().SystIndex(ShortName());

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
    u.i0 = uo.ClosestShiftIndex(ShortName(), x, ESide::kBelow, &x0);
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
  /*
  const std::vector<const ISyst*>& GetSBNGenieWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(ret.empty()){
      const UniverseOracle& uo = UniverseOracle::Instance();
      ret.reserve(uo.Systs().size());
      for(const std::string& name: uo.Systs()){
        if(name.find("genie") == std::string::npos) continue;

        if(name.find("NonRes") == 0){
          // These NonRes parameters need special handling. There are two
          // copies, supposed to apply to CC and NC respectively, but instead
          // the weights are filled for all events.
          if(name.find("Alt") == std::string::npos){
            ret.push_back(new SBNWeightSyst(name, kIsNC));
          }
          else{
            ret.push_back(new SBNWeightSyst(name, kIsCC));
          }
        }
        else{
          // Regular case
          ret.push_back(new SBNWeightSyst(name));
        }
      }
    }
    return ret;
  }
  */

  // --------------------------------------------------------------------------
  /*
  const std::vector<const ISyst*>& GetSBNFluxWeightSysts()
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
      const std::vector<const ISyst*>& f = GetSBNFluxWeightSysts();
      ret.reserve(g.size()+f.size());
      ret.insert(ret.end(), g.begin(), g.end());
      ret.insert(ret.end(), f.begin(), f.end());
    }
    return ret;
  }
  */
}
