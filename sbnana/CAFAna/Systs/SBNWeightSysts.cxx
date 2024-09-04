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
  double UniverseWeight::operator()(const caf::SRTrueInteractionProxy* nu) const
  {

    if(fPSetIdx == -1){
      const UniverseOracle& uo = UniverseOracle::Instance();
      fPSetIdx = uo.ParameterSetIndex(fPSetName);
    }

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = nu->wgt;
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
    this->Shift(x, &sr->truth, weight);
  }

  // --------------------------------------------------------------------------
  void SBNWeightSyst::Shift(double x, caf::SRTrueInteractionProxy* nu, double& weight) const
  {
    if(nu->index < 0) return;

    if(fIdx == -1) fIdx = UniverseOracle::Instance().SystIndex(ShortName());

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = nu->wgt;
    if(wgts.empty()) return;

    // Check if 0-to-1 dial
    // For 0-to-1 dials, e.g., GENIE DecayAngMEC,
    // we only save one universe with x=1.
    // In this case, 
    bool IsMorphDial = false;
    const std::vector<float>& v = UniverseOracle::Instance().ShiftsForSyst( ShortName() );
    if(v.size()==1){
      if( fabs(v[0] - 1.0) < 1E-5 ){
        IsMorphDial = true;
      }
    }

    if(IsMorphDial){
      double fullrw = wgts[fIdx].univ[0];
      double this_rw = x * fullrw + (1.-x) * 1.;
      weight *= this_rw;
    }
    else{
      const Univs u = GetUnivs(x);

      double y = 0;
      if(u.w0 != 0) y += u.w0 * wgts[fIdx].univ[u.i0];
      if(u.w1 != 0) y += u.w1 * wgts[fIdx].univ[u.i1];

      weight *= y;
    }

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
    return {"AhtBY",
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
    return {"expskin",
            "horncurrent",
            "nucleoninexsec",
            "nucleonqexsec",
            "nucleontotxsec",
            "pioninexsec",
            "pionqexsec",
            "piontotxsec"};
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
  
  // --------------------------------------------------------------------------
  SBNWeightMirrorSyst::SBNWeightMirrorSyst(const std::string& systName)
    : SBNWeightSyst(systName)
  {
  }

  // --------------------------------------------------------------------------
  void SBNWeightMirrorSyst::Shift(double x, caf::SRSliceProxy* sr, double& weight) const
  {
    this->Shift(x, &sr->truth, weight);
  }

  // --------------------------------------------------------------------------
  void SBNWeightMirrorSyst::Shift(double x, caf::SRTrueInteractionProxy* nu, double& weight) const
  {

    double mirrored_x = x<0 ? -x : x;

    if(nu->index < 0) return;

    if(fIdx == -1) fIdx = UniverseOracle::Instance().SystIndex(ShortName());

    const caf::Proxy<std::vector<caf::SRMultiverse>>& wgts = nu->wgt;
    if(wgts.empty()) return;

    // Find dial=1sigma
    const std::vector<float>& v = UniverseOracle::Instance().ShiftsForSyst( ShortName() );
    int iu_onesig = -1;
    for(unsigned int iu=0; iu<v.size(); iu++){
      if( fabs(v[iu] - 1.0) < 1E-5 ){
        iu_onesig = iu;
        break;
      }
    }

    if(iu_onesig>=0){
      double fullrw = wgts[fIdx].univ[iu_onesig];
      double this_rw = mirrored_x * fullrw + (1.-mirrored_x) * 1.;
      weight *= this_rw;
    }
    else{
      const Univs u = GetUnivs(mirrored_x);

      double y = 0;
      if(u.w0 != 0) y += u.w0 * wgts[fIdx].univ[u.i0];
      if(u.w1 != 0) y += u.w1 * wgts[fIdx].univ[u.i1];

      weight *= y;

    }
  }

  // --------------------------------------------------------------------------
  SBNWeightAbsVarSyst::SBNWeightAbsVarSyst(const std::string& systName, std::vector< std::pair<double, double> > _absvar_to_sigma)
    : SBNWeightSyst(systName), absvar_to_sigma(_absvar_to_sigma)
  {

    assert( absvar_to_sigma.size()>=2 );

    // Sort the vector by the sigma values (second in the pair)
    // first: absolute value, second: sigma
    std::sort(
      absvar_to_sigma.begin(),
      absvar_to_sigma.end(),
      [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.second < b.second;
      }
    );

  }

  double SBNWeightAbsVarSyst::ConvertToAbsolute(double x) const {

    // x: a given sigma
    // Use absvar_to_sigma to decide absvar

    static size_t n = absvar_to_sigma.size();

    // absvar_to_sigma
    // - first: absolute value
    // - second: sigma

    // index for the absvar_to_sigma vector
    int idx_left = 0;
    int idx_right = 1;

    // 1) When the given sigma x is lower than the left-boundary,
    //    Use the first two values, 0 and 1, and extrapolate
    if (x <= absvar_to_sigma[0].second){
      idx_left = 0;
      idx_right = 1;
    }
    // 2) When the given sigma x is higher than the right-boundary,
    //    Use the last two values, n-2 and n-1, and extrapolate
    else if(x >= absvar_to_sigma[n-1].second) {
      idx_left = n-2;
      idx_right = n-1;
    }
    // 3) if in range, find the closest neighbor
    else{

      for (size_t i = 0; i < absvar_to_sigma.size() - 1; ++i) {
        if (x >= absvar_to_sigma[i].second && x <= absvar_to_sigma[i+1].second) {
          idx_left = i;
          idx_right = i+1;
          break;
        }
      }

    }

    return absvar_to_sigma[idx_left].first +
           (x - absvar_to_sigma[idx_left].second) *
           (absvar_to_sigma[idx_right].first - absvar_to_sigma[idx_left].first) /
           (absvar_to_sigma[idx_right].second - absvar_to_sigma[idx_left].second);

  }

  // --------------------------------------------------------------------------
  void SBNWeightAbsVarSyst::Shift(double x, caf::SRSliceProxy* sr, double& weight) const
  {
    this->Shift(x, &sr->truth, weight);
  }

  // --------------------------------------------------------------------------
  void SBNWeightAbsVarSyst::Shift(double x, caf::SRTrueInteractionProxy* nu, double& weight) const
  {

    if(nu->index < 0) return;

    // convert x (sigma) into an absoluate value which is saved as univ[] in CAF
    double this_abs = ConvertToAbsolute(x);

    return SBNWeightSyst::Shift(this_abs, nu, weight);

  }

  // --------------------------------------------------------------------------
  SBNWeightUnivSyst::SBNWeightUnivSyst(
    const std::string& systName,
    const std::string& psetName,
    std::vector< std::pair<unsigned, double> > _univ_to_sigma
  ) : SBNWeightSyst(systName), kPsetName(psetName), univ_to_sigma(_univ_to_sigma)
  {

    assert( univ_to_sigma.size()>=2 );

    // Sort the vector by the sigma values (second in the pair)
    // first: univ, second: sigma
    std::sort(
      univ_to_sigma.begin(),
      univ_to_sigma.end(),
      [](const std::pair<unsigned, double>& a, const std::pair<unsigned, double>& b) {
        return a.second < b.second;
      }
    );

    for(const auto& it: univ_to_sigma){
      std::cout << "[SBNWeightUnivSyst::SBNWeightUnivSyst] GetTruthUniverseWeight with " << kPsetName << ", " << it.first << std::endl;
      vec_TruthVar.push_back( GetTruthUniverseWeight(kPsetName, it.first) );
    }

  }

  // --------------------------------------------------------------------------
  void SBNWeightUnivSyst::Shift(double x, caf::SRSliceProxy* sr, double& weight) const
  {
    this->Shift(x, &sr->truth, weight);
  }

  // --------------------------------------------------------------------------
  void SBNWeightUnivSyst::Shift(double x, caf::SRTrueInteractionProxy* nu, double& weight) const
  {

    if(nu->index < 0) return;

    // x: a given sigma
    // Use absvar_to_sigma to decide absvar

    static size_t n = univ_to_sigma.size();

    // absvar_to_sigma
    // - first: absolute value
    // - second: sigma

    int idx_left = 0;
    int idx_right = 1;

    // 1) When the given sigma x is lower than the left-boundary,
    //    Use the first two values, 0 and 1, and extrapolate
    if (x <= univ_to_sigma[0].second){
      idx_left = 0;
      idx_right = 1;
    }
    // 2) When the given sigma x is higher than the right-boundary,
    //    Use the last two values, n-2 and n-1, and extrapolate
    else if(x >= univ_to_sigma.back().second) { 
      idx_left = n-2;
      idx_right = n-1;
    }
    // 3) if in range, find the closest neighbor
    else{

      for (size_t i = 0; i < univ_to_sigma.size() - 1; ++i) { 
        if (x >= univ_to_sigma[i].second && x <= univ_to_sigma[i+1].second) {
          idx_left = i;
          idx_right = i+1;
          break;
        }
      }

    }

    // map_univ_to_TruthVar
    double rw_left = vec_TruthVar[idx_left](nu);
    double rw_right = vec_TruthVar[idx_right](nu);
    double sigma_left = univ_to_sigma[idx_left].second;
    double sigma_right = univ_to_sigma[idx_right].second;
    double this_rw = rw_left + (x-sigma_left) * (rw_left - rw_right) / (sigma_left - sigma_right);
    if(this_rw<0) this_rw = 0.;
    weight *= this_rw;

  }

}
