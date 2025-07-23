#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "cetlib/search_path.h"

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

const Var kIcarus202401RecoENu = kIcarus202401RecoMuonE + kIcarus202401RecoProtonKE + kIcarus202401RecoPionE;

static bool Icarus202401_contained(const caf::SRTrackProxy& trk)
{
  return ((trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5)
             || (trk.end.x > 61.94 + 5 && trk.end.x < 358.49 - 5))
           && trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5
           && trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5;
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


const Var kMuon_length_reco([](const caf::SRSliceProxy* slc){
        double E =0;
 int mu_idx = kIcarus202401MuonIdx(slc);
  if(mu_idx < 0) return 0.;
  const auto &muTrk = slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk;
  if(Icarus202401_contained(muTrk)==false)return 0.;
    if(mu_idx != -1){E = slc->reco.pfp[mu_idx].trk.len;}
    return E;
});


const Var kProton_length_reco([](const caf::SRSliceProxy* slc){
        //int ipfp_mu = -1;
        //int ipfp_pro = -1;
        double E =0;
int mu_idx = kIcarus202401MuonIdx(slc);

        double max_length = 0;
        for ( std::size_t ipfp(0); ipfp < slc->reco.npfp ; ++ipfp ){
            if(int(ipfp)==mu_idx)continue;
             auto const& trk = slc->reco.pfp[ipfp].trk;
            const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
            if(Icarus202401_proton_cut(slc->reco.pfp[ipfp].trk) && slc->reco.pfp[ipfp].parent_is_primary && Atslc < 10.0 && std::hypot(slc->reco.pfp[ipfp].trk.rangeP.p_proton, mp) - mp > 0.05){
                if(max_length<slc->reco.pfp[ipfp].trk.len){max_length=slc->reco.pfp[ipfp].trk.len;}
            }
                    }

        E = max_length;
    return E;
});

}






