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

        // These ones have to be handled specially
        const std::string prefix = "NonResR";

        if(name.find(prefix) == 0){
          // This logic is copied from GetSBNGenieWeightSysts()
          for(std::string npi: {"1pi", "2pi"}){
            // They need to be split into CC and NC variants
            for(std::string ccncStr: {"CC", "NC"}){
              const Cut& ccncCut = (ccncStr == "CC") ? kIsCC : kIsNC;

              // And restricted to either neutrino or antineutrino, taking
              // advantage of symmetries to implement both proton and neutron
              // versions.

              // Originals
              if(name == prefix + "vp" + npi + ccncStr){
                fSystIdxs.push_back(uo.SystIndex(prefix + "vp" + npi));
                fUnivOffsets.push_back(0);
                fUnivCuts.push_back(ccncCut && !kIsAntiNu);
              }
              if(name == prefix + "vbarp" + npi + ccncStr){
                fSystIdxs.push_back(uo.SystIndex(prefix + "vbarp" + npi));
                fUnivOffsets.push_back(1);
                fUnivCuts.push_back(ccncCut &&  kIsAntiNu);
              }

              // Switched variants
              if(name == prefix + "vbarn" + npi + ccncStr){
                fSystIdxs.push_back(uo.SystIndex(prefix + "vp" + npi));
                fUnivOffsets.push_back(2);
                fUnivCuts.push_back(ccncCut &&  kIsAntiNu);
              }
              if(name == prefix + "vn" + npi + ccncStr){
                fSystIdxs.push_back(uo.SystIndex(prefix + "vbarp" + npi));
                fUnivOffsets.push_back(3);
                fUnivCuts.push_back(ccncCut && !kIsAntiNu);
              }
            }
          }
        }
        else{
          // "Normal" syst
          fSystIdxs.push_back(uo.SystIndex(name));
          fUnivOffsets.push_back(0);
          fUnivCuts.push_back(kNoCut);
        }
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
                               const std::string& knobName,
                               const SliceCut& cut)
    : ISyst(systName, systName),
      fIdx(-1),
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
  const std::vector<const ISyst*>& GetSBNGenieWeightSysts()
  {
    static std::vector<const ISyst*> ret;
    if(!ret.empty()) return ret;

    // We can't ask the UniverseOracle about this, because it doesn't get
    // properly configured until it's seen its first CAF file.
    const std::vector<std::string> names = {"DISAth",
                                            "DISBth",
                                            "DISCv1u",
                                            "DISCv2u",
                                            "IntraNukeNabs",
                                            "IntraNukeNcex",
                                            "IntraNukeNinel",
                                            "IntraNukeNmfp",
                                            "IntraNukeNpi",
                                            "IntraNukePIabs",
                                            "IntraNukePIcex",
                                            "IntraNukePIinel",
                                            "IntraNukePImfp",
                                            "IntraNukePIpi",
                                            "NC",
                                            "ResDecayGamma",
                                            "CCResAxial",
                                            "CCResVector",
                                            // Disabled due to NaNs
                                            // "CohMA",
                                            // "CohR0",
                                            "NCELaxial",
                                            "NCELeta",
                                            "NCResAxial",
                                            "NCResVector",
                                            "QEMA"};

    for(const std::string& name: names) ret.push_back(new SBNWeightSyst(name));

    // These ones have to be handled specially
    const std::string prefix = "NonResR";

    for(std::string npi: {"1pi", "2pi"}){
      // They need to be split into CC and NC variants
      for(std::string ccncStr: {"CC", "NC"}){
        const Cut& ccncCut = (ccncStr == "CC") ? kIsCC : kIsNC;

        // And restricted to either neutrino or antineutrino, taking advantage
        // of symmetries to implement both proton and neutron versions.

        // First constructor argument is the syst name, second the knob to use

        // Originals
        ret.push_back(new SBNWeightSyst(prefix + "vp"    + npi + ccncStr,
                                        prefix + "vp"    + npi,
                                        ccncCut && !kIsAntiNu));
        ret.push_back(new SBNWeightSyst(prefix + "vbarp" + npi + ccncStr,
                                        prefix + "vbarp" + npi,
                                        ccncCut &&  kIsAntiNu));

        // Switched variants
        ret.push_back(new SBNWeightSyst(prefix + "vbarn" + npi + ccncStr,
                                        prefix + "vp"    + npi,
                                        ccncCut &&  kIsAntiNu));
        ret.push_back(new SBNWeightSyst(prefix + "vn"    + npi + ccncStr,
                                        prefix + "vbarp" + npi,
                                        ccncCut && !kIsAntiNu));
      }
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
