#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "cetlib/search_path.h"

#include "TVector3.h"
#include "Math/Vector3D.h"
#include "Math/VectorUtil.h"

namespace ana {

namespace chi2pid {

void Chi2PID::check_graphs() {
  if(dedx_range_pro == NULL || dedx_range_ka == NULL ||
     dedx_range_pi == NULL || dedx_range_mu == NULL) {
    cet::search_path sp("FW_SEARCH_PATH");

    std::string kdEdXUncTemplateFileName = "template_dEdXUncertainty.root";
    std::string kdEdXUncTemplateFullFilePath;

    if (!sp.find_file(kdEdXUncTemplateFileName, kdEdXUncTemplateFullFilePath)) {
      std::cerr << "ICARUS 202401 numu selection: failed to locate dE/dx template file '"
        << kdEdXUncTemplateFileName << "'" << std::endl;
      std::abort();
    }

    file_dEdXUncTemplate = TFile::Open(kdEdXUncTemplateFullFilePath.c_str());
    dedx_unc_template = file_dEdXUncTemplate->Get<TGraph2D>("dEdXRelUncertainty_dEdX_vs_phi");

    std::string kChi2TemplateFileName = "dEdxrestemplates.root";
    std::string kChi2TemplateFullFilePath;

    if (!sp.find_file(kChi2TemplateFileName, kChi2TemplateFullFilePath)) {
      std::cerr << "\nICARUS 202401 numu selection: failed to locate dE/dx template file '"
        << kChi2TemplateFileName << "'" << std::endl;
      std::abort();
    }

    file_Chi2Template = TFile::Open(kChi2TemplateFullFilePath.c_str());
    dedx_range_pro = file_Chi2Template->Get<TProfile>("dedx_range_pro");
    dedx_range_ka  = file_Chi2Template->Get<TProfile>("dedx_range_ka");
    dedx_range_pi  = file_Chi2Template->Get<TProfile>("dedx_range_pi");
    dedx_range_mu  = file_Chi2Template->Get<TProfile>("dedx_range_mu");

  }
}

Chi2PID chi2_calculator;

}

static constexpr float mmu = 0.106f;
static constexpr float mp  = 0.9383f;
static constexpr float mpi = 0.13957f;

// Mostly just copied from NumuVarsIcarus202106.cxx
const Var kIcarus202401MuonIdx([](const caf::SRSliceProxy* slc) -> int {
      // The (dis)qualification of a slice is based upon the track level information.
      float Longest(0);
      int PTrackInd(-1);
      for (std::size_t i(0); i < slc->reco.npfp; ++i) {
        auto const& pfp = slc->reco.pfp.at(i);
        if (pfp.trackScore < 0.5) { continue; }
	//if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x)){ continue; }
        auto const& trk = pfp.trk;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

        //int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
        int plane = 2; // Hard code collection plane for now since induction 2 has peak at higher chi2

        //float Chi2Proton = trk.chi2pid[plane].chi2_proton;
        //float Chi2Muon = trk.chi2pid[plane].chi2_muon;
        auto chi2 = chi2pid::chi2_calculator.calculate_chi2(trk.calo[plane]);
        float Chi2Proton = chi2.chi2_proton;
        float Chi2Muon = chi2.chi2_muon;
  
        const bool Contained = ( !isnan(trk.end.x) &&
            (( trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5 ) ||
             ( trk.end.x > 61.94 + 5 && trk.end.x < +358.49 - 5 )) &&
            !isnan(trk.end.y) &&
            ( trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5 ) &&
            !isnan(trk.end.z) &&
            ( trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5 ) );
        const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest ) {
          Longest = trk.len;
          PTrackInd = i;
        }
      }
      return PTrackInd;
});

static bool Icarus202401_proton_cut(const caf::SRTrackProxy& trk)
{
  //int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
  int plane = 2; // Hard code collection plane for now since induction 2 has peak at higher chi2
  auto chi2 = chi2pid::chi2_calculator.calculate_chi2(trk.calo[plane]);
  return chi2.chi2_proton < 100 && chi2.chi2_proton > 0;
}

static bool Icarus202401_proton_cut_exist(const caf::SRTrackProxy& trk)
{ 
  //int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
    int plane = 2; // Hard code collection plane for now since induction 2 has peak at higher chi2
    auto chi2 = chi2pid::chi2_calculator.calculate_chi2(trk.calo[plane]);
    return chi2.chi2_proton > 0;
}

const Var kIcarus202401NumPions([](const caf::SRSliceProxy* slc)
{
  int count = 0;
  auto idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    if (pfp.trackScore < 0.5) { continue; }
    auto const& trk = pfp.trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
    //if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
    if(pfp.id != muID && !Icarus202401_proton_cut(trk) && Icarus202401_proton_cut_exist(trk) && AtSlice && std::hypot(pfp.trk.rangeP.p_pion, mpi) - mpi > 0.025)
      ++count;
  }
  return count;
});

template<int pdg, int threshold>
const MultiVar kTruePrimaryPs([](const caf::SRSliceProxy* slc){
    std::vector<double> ps;
    for(size_t i = 0; i < slc->truth.prim.size(); i++) {
        const auto& p = slc->truth.prim[i];
        double visE = 0;
        if(!std::isnan(p.plane[0][2].visE)) visE += p.plane[0][2].visE;
        if(!std::isnan(p.plane[1][2].visE)) visE += p.plane[1][2].visE;
        if(std::abs(p.pdg) == pdg && visE >= threshold/1000.0) ps.push_back(i);
    }
    return ps;
});

const Var kIcarus202401NumProtons([](const caf::SRSliceProxy* slc)
{
  int count = 0;
  auto idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    if (pfp.trackScore < 0.4) { continue; }
    //if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
    auto const& trk = pfp.trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
    if(pfp.id != muID && Icarus202401_proton_cut(trk) && AtSlice && std::hypot(pfp.trk.rangeP.p_proton, mp) - mp > 0.05)
      ++count;
  }
  return count;
});

const Var kIcarus202401NumShowers([](const caf::SRSliceProxy* slc){
  int count = 0;
  for(const auto& pfp: slc->reco.pfp) {
    [[maybe_unused]]
    const float Atslc = std::hypot(slc->vertex.x - pfp.shw.start.x,
                                   slc->vertex.y - pfp.shw.start.y,
                                   slc->vertex.z - pfp.shw.start.z);
    const bool AtSlice = ( Atslc < 50.0 && pfp.parent_is_primary);
    const bool isShower = pfp.trackScore > 0 && (pfp.trackScore < 0.4 || (pfp.trackScore < 0.5 && !Icarus202401_proton_cut(pfp.trk)));
    if(isShower && AtSlice && pfp.shw.plane[2].energy > 0.025) count++;
  }
  return count;
});

const Var kIcarus202401RecoMuonE([](const caf::SRSliceProxy* slc){
  return std::hypot(kIcarus202401RecoMuonP(slc), mmu);
});

const Var kIcarus202401RecoProtonKE([](const caf::SRSliceProxy* slc){
  double E = 0;
  for(const auto P: kIcarus202401RecoProtonP(slc))
    E += std::hypot(P, mp) - mp;
  return E;
});

const Var kIcarus202401RecoPionE([](const caf::SRSliceProxy* slc){
  double E = 0;
  for(const auto P: kIcarus202401RecoPionP(slc))
    E += std::hypot(P, mpi);
  return E;
});

const Var kIcarus202401RecoBindingE([](const caf::SRSliceProxy *slc) {
  return 0.0309 * kIcarus202401RecoProtonP(slc).size();
});

const Var kIcarus202401RecoENu = kIcarus202401RecoMuonE + kIcarus202401RecoProtonKE + kIcarus202401RecoPionE + kIcarus202401RecoBindingE;

static bool Icarus202401_contained(const caf::SRTrackProxy& trk)
{
  const double x = trk.end.x, y = trk.end.y, z = trk.end.z;
  if(std::isnan(x) || std::isnan(y) || std::isnan(z)) return false;

  bool x_contained_E = x <  -61.94 - 5 && x > -358.49 + 5;
  bool x_contained_W = x >   61.94 + 5 && x <  358.49 - 5;
  bool y_contained   = y > -181.86 + 5 && y <  134.96 - 5;
  bool z_contained   = z > -894.95 + 5 && z <  894.95 - 5;

  // Check triangles at corners of the detector
  bool in_corner = y <  1.732007 * z - 1687.5114 
                || y > -1.732007 * z + 1640.6114
                || y >  1.732007 * z + 1640.6114 
                || y < -1.732007 * z - 1687.5114;

  return (x_contained_E || x_contained_W) && y_contained && z_contained 
         && !in_corner;
}

const Var kIcarus202401RecoMuonP([](const caf::SRSliceProxy *slc){
  double ret = 0;
  if(kIcarus202401MuonIdx(slc) < 0) return ret;
  const auto &muTrk = slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk;
  ret = Icarus202401_contained(muTrk) ? muTrk.rangeP.p_muon : 
               muTrk.mcsP.is_bwd_muon ? muTrk.mcsP.bwdP_muon : 
                                        muTrk.mcsP.fwdP_muon;
  return ret;
});

const MultiVar kIcarus202401RecoProtonP([](const caf::SRSliceProxy* slc){
  std::vector<double> Ps;
  int mu_idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if(mu_idx >= 0) muID = slc->reco.pfp[mu_idx].id;
  for(const auto& pfp: slc->reco.pfp) {
    if(pfp.id == muID) continue;
    if(pfp.trackScore < 0.4) continue;
    auto const& trk = pfp.trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    if(Icarus202401_proton_cut(pfp.trk) && pfp.parent_is_primary && Atslc < 10.0 && std::hypot(pfp.trk.rangeP.p_proton, mp) - mp > 0.05)
      Ps.push_back(pfp.trk.rangeP.p_proton);
  }
  return Ps;
});

const MultiVar kIcarus202401RecoPionP([](const caf::SRSliceProxy* slc){
  std::vector<double> Ps;
  int mu_idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if(mu_idx >= 0) muID = slc->reco.pfp[mu_idx].id;
  for(const auto& pfp: slc->reco.pfp) {
    if(pfp.id == muID) continue;
    if(pfp.trackScore < 0.5) continue;
    auto const& trk = pfp.trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    if(Atslc < 10.0 && !Icarus202401_proton_cut(pfp.trk) && pfp.parent_is_primary && Icarus202401_proton_cut_exist(trk) && std::hypot(pfp.trk.rangeP.p_pion, mpi) - mpi > 0.025)
      Ps.push_back(pfp.trk.rangeP.p_pion);
  }
  return Ps;
});

const Var kIcarus202401RecoTransP([](const caf::SRSliceProxy *slc){
  double ret = 0;
  double p_mu_x = 0; double p_p_x = 0;
  double p_mu_y = 0; double p_p_y = 0;
  int mu_idx = kIcarus202401MuonIdx(slc);
  if(mu_idx < 0) return p_mu_x;
  const auto &muTrk = slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk;
  if(Icarus202401_contained(muTrk)==false)return p_mu_x;
  p_mu_x =muTrk.rangeP.p_muon*muTrk.dir.x;
  p_mu_y =muTrk.rangeP.p_muon*muTrk.dir.y;
  //float mp = 0.9383;
  int muID = -1;
  if(mu_idx >= 0) muID = slc->reco.pfp[mu_idx].id;
  for(const auto& pfp: slc->reco.pfp) {
    if(pfp.id == muID) continue;
    if(pfp.trackScore < 0.4) continue;
    auto const& trk = pfp.trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    if(Icarus202401_proton_cut(pfp.trk) && pfp.parent_is_primary && Atslc < 10.0 && std::hypot(pfp.trk.rangeP.p_proton, mp) - mp > 0.05){
        p_p_x+=pfp.trk.rangeP.p_proton*pfp.trk.dir.x;
        p_p_y+=pfp.trk.rangeP.p_proton*pfp.trk.dir.y;
    }
  }

  double p_tot_x = p_mu_x + p_p_x;
  double p_tot_y = p_mu_y + p_p_y;
  ret = std::hypot(p_tot_x, p_tot_y);

  return ret;
});

const Var kICARUS202401TrueLeadingProtonMomentum([] (const caf::SRSliceProxy* slc) {
  double protonE = -9999.;
  double protonMom = -9999.;
  for (auto& prim: slc->truth.prim) {
    if (prim.pdg==2212) {
      if (prim.startE > protonE) {
        protonMom = sqrt(prim.startp.x*prim.startp.x + prim.startp.y*prim.startp.y + prim.startp.z*prim.startp.z);
        protonE = prim.startE;
      }
    }
  }
  return protonMom;
});

const Var kICARUS202401TrueLeadingProtonCosTheta([] (const caf::SRSliceProxy* slc) {
  double protonE = -9999.;
  int nproton = 0;
  //Get neutrino direction
  ROOT::Math::XYZVector nuvec = ROOT::Math::XYZVector(slc->truth.momentum.x,slc->truth.momentum.y,slc->truth.momentum.z);
  ROOT::Math::XYZVector pvec(0.,0.,0.);
  for (auto& prim: slc->truth.prim) {
    if (prim.pdg==2212) {
      nproton++;
      if (prim.startE > protonE) {
	protonE = prim.startE;
	pvec = ROOT::Math::XYZVector(prim.startp.x,prim.startp.y,prim.startp.z);
      }
    }
  }
  double protoncosth = ROOT::Math::VectorUtil::CosTheta(nuvec,pvec);
  if (nproton < 1) {
    protoncosth = -9999.;
  }
  return protoncosth;
});

const Var kIcarus202401TrueMuonP([](const caf::SRSliceProxy *slc) -> double {
  for(const auto &p: slc->truth.prim)
    if(std::abs(p.pdg) == 13)
      return std::hypot(p.startp.x, p.startp.y, p.startp.z);
  return -1.0;
});

const Var kIcarus202401TrueMuonCosTheta([](const caf::SRSliceProxy *slc) -> double {
  for(const auto &p: slc->truth.prim)
    if(std::abs(p.pdg) == 13)
      return p.startp.z/std::hypot(p.startp.x, p.startp.y, p.startp.z);
  return -2.0;
});

static double Icarus202401_T3D_angle_mup ( const caf::SRSliceProxy* islc, int ipfp_mu, int ipfp_pro ) { 

    float p_mu_x,p_mu_y,p_mu_z;
    float p_p_x,p_p_y,p_p_z;
    //float p_tot_x,p_tot_y,p_tot_z;

    //Muon momentum
    p_mu_x=(islc->reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc->reco.pfp[ipfp_mu].trk.dir.x; //GeV
    p_mu_y=(islc->reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc->reco.pfp[ipfp_mu].trk.dir.y;
    p_mu_z=(islc->reco.pfp[ipfp_mu].trk.rangeP.p_muon)*islc->reco.pfp[ipfp_mu].trk.dir.z;
    TVector3 mu_vector(p_mu_x,p_mu_y,p_mu_z);

    //Proton momentum
    p_p_x =(islc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc->reco.pfp[ipfp_pro].trk.dir.x;
    p_p_y =(islc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc->reco.pfp[ipfp_pro].trk.dir.y;
    p_p_z =(islc->reco.pfp[ipfp_pro].trk.rangeP.p_proton)*islc->reco.pfp[ipfp_pro].trk.dir.z;
    TVector3 pro_vector(p_p_x,p_p_y,p_p_z);

    return cos(mu_vector.Angle(pro_vector)); 
}


const Var kIcarus202401RecoLeadingProtonLen([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -1.0;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.len;
    return -1.0;
});

const Var kIcarus202401RecoLeadingProtonP([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return 0.0;
    return *std::max_element(ps.begin(), ps.end());
});

const Var kIcarus202401RecoMuonLen([](const caf::SRSliceProxy* islc)->double{
    int ipfp_mu = -1;
    ipfp_mu =kIcarus202401MuonIdx(islc);
    return ipfp_mu != -1 ? (double)islc->reco.pfp[ipfp_mu].trk.len : -1000.0;
});

const Var kIcarus202401RecoMuonDirX([](const caf::SRSliceProxy* islc)->double{
    int ipfp_mu = -1;
    ipfp_mu =kIcarus202401MuonIdx(islc);
    return ipfp_mu != -1 ? (double)islc->reco.pfp[ipfp_mu].trk.dir.x : -1000.0;
});

const Var kIcarus202401RecoMuonDirY([](const caf::SRSliceProxy* islc)->double{
    int ipfp_mu = -1;
    ipfp_mu =kIcarus202401MuonIdx(islc);
    return ipfp_mu != -1 ? (double)islc->reco.pfp[ipfp_mu].trk.dir.y : -1000.0;
});

const Var kIcarus202401RecoMuonDirZ([](const caf::SRSliceProxy* islc)->double{
    int ipfp_mu = -1;
    ipfp_mu =kIcarus202401MuonIdx(islc);
    return ipfp_mu != -1 ? (double)islc->reco.pfp[ipfp_mu].trk.dir.z : -1000.0;
});

const Var kIcarus202401RecoLeadingProtonDirX([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -2.0;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.dir.x;
    return -2.0;
});

const Var kIcarus202401RecoLeadingProtonDirY([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -2.0;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.dir.y;
    return -2.0;
});

const Var kIcarus202401RecoLeadingProtonDirZ([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -2.0;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.dir.z;
    return -2.0;
});

const Var kIcarus202401RecoLeadingProtonCosTheta([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -2.0;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.costh;
    return -2.0;
});

const Var kIcarus202401RecoLeadingProtonPhi([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -5.0;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.phi;
    return -5.0;
});

const Var kIcarus202401RecoLeadingProtonEndX([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -100000;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.end.x;
    return -100000;;
});

const Var kIcarus202401RecoLeadingProtonEndY([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -100000;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.end.x;
    return -100000;;
});

const Var kIcarus202401RecoLeadingProtonEndZ([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -100000;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return pfp.trk.end.x;
    return -100000;;
});

const Var kIcarus202401RecoMuLeadPOpeningAngle([](const caf::SRSliceProxy* islc)->double{
        int ipfp_mu = -1;
        double E =0;
        ipfp_mu =kIcarus202401MuonIdx(islc);
        const auto ps = kIcarus202401RecoProtonP(islc);
        if(!ps.size()) return -1000.0;
        double p = *std::max_element(ps.begin(), ps.end());
        size_t ipfp_pro = 0;
        for(; ipfp_pro < islc->reco.pfp.size(); ++ipfp_pro)
          if(islc->reco.pfp[ipfp_pro].trk.rangeP.p_proton == p)
            break;
        E = Icarus202401_T3D_angle_mup(islc,ipfp_mu,ipfp_pro);
    return E;
});

const Var kIcarus202401MuonChi2Mu([](const caf::SRSliceProxy* islc) ->double {
        int ipfp_mu = -1;
        ipfp_mu =kIcarus202401MuonIdx(islc);
        if(ipfp_mu < 0) return -1;
        auto chi2 = chi2pid::chi2_calculator.calculate_chi2(islc->reco.pfp[ipfp_mu].trk.calo[2]);
        return chi2.chi2_muon;
});

const Var kIcarus202401LeadingProtonChi2Proton([](const caf::SRSliceProxy* islc)->double{
    const auto ps = kIcarus202401RecoProtonP(islc);
    if(!ps.size()) return -1.0;
    double p = *std::max_element(ps.begin(), ps.end());
    for(const auto &pfp: islc->reco.pfp)
      if(pfp.trk.rangeP.p_proton == p)
        return chi2pid::chi2_calculator.calculate_chi2(pfp.trk.calo[2]).chi2_proton;
    return -1.0;
});

const Var kIcarus202401RecoMuonCosTheta([](const caf::SRSliceProxy *slc)->double{
  if(kIcarus202401MuonIdx(slc) < 0) return -2.0;
  return slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk.costh;
});
const Var kIcarus202401RecoMuonPhi([](const caf::SRSliceProxy *slc)->double{
  if(kIcarus202401MuonIdx(slc) < 0) return -2.0;
  return slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk.phi;
});
const Var kIcarus202401RecoMuonEndX([](const caf::SRSliceProxy *slc)->double{
  if(kIcarus202401MuonIdx(slc) < 0) return -9999.;
  return slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk.end.x;
});
const Var kIcarus202401RecoMuonEndY([](const caf::SRSliceProxy *slc)->double{
  if(kIcarus202401MuonIdx(slc) < 0) return -9999.;
  return slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk.end.y;
});
const Var kIcarus202401RecoMuonEndZ([](const caf::SRSliceProxy *slc)->double{
  if(kIcarus202401MuonIdx(slc) < 0) return -9999.;
  return slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk.end.z;
});
const Var kIcarus202401BaryFMDeltaZ = SIMPLEVAR(barycenterFM.deltaZ_Trigger);

}






