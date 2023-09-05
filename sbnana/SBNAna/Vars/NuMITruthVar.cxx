#include "sbnana/SBNAna/Vars/NuMITruthVar.h"

namespace ana{

  // Neutrino/interaction

  const TruthVar kTruth_NProton_Primary([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 2212 ) NPtl++;
    }
    return NPtl;
  });
  const TruthVar kTruth_NNeutron_Primary([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 2112 ) NPtl++;
    }
    return NPtl;
  });
  const TruthVar kTruth_Npip_Primary([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 211 ) NPtl++;
    }
    return NPtl;
  });
  const TruthVar kTruth_Npip_All([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( prim.pdg == 211 ) NPtl++;
    }
    return NPtl;
  });
  const TruthVar kTruth_Npim_Primary([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == -211 ) NPtl++;
    }
    return NPtl;
  });
  const TruthVar kTruth_Npim_All([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( prim.pdg == -211 ) NPtl++;
    }
    return NPtl;
  });
  const TruthVar kTruth_Npi0_Primary([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( abs(prim.pdg) == 111 ) NPtl++;
    }
    return NPtl;
  });
  const TruthVar kTruth_Npi0_All([](const caf::SRTrueInteractionProxy *nu) -> int {
    int NPtl = 0;
    for ( auto const& prim : nu->prim ) {
      if ( abs(prim.pdg) == 111 ) NPtl++;
    }
    return NPtl;
  });

  // Muon

  const TruthVar kTruth_MuonIndex([](const caf::SRTrueInteractionProxy *nu) -> int {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < nu->prim.size(); ++i){
      // primary
      if( nu->prim.at(i).start_process!=0 ) continue;
      // muon
      if( abs(nu->prim.at(i).pdg)!=13 ) continue;
      // non-nan genE
      if(isnan(nu->prim.at(i).genE)) continue;

      double this_E = nu->prim.at(i).genE;
      // if larger E, update
      if(this_E>max_E){
        max_E = this_E;
        truth_idx = i;
      }
    }
    return truth_idx;
  });
  const TruthVar kTruth_MuonNuCosineTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_MuonIndex(nu);
    if(truth_idx>=0){

      const auto& p_mu = nu->prim.at(truth_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      ret = vec_p_mu.Unit().Dot( vec_p_nu.Unit() );

    }

    return ret;
  });
  const TruthVar kTruth_MuonCosThBeam([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_MuonIndex(nu);
    if(truth_idx>=0){

      TVector3 rFromNuMI(315.120380, 33.644912, 733.632532);

      const auto& p_mu = nu->prim.at(truth_idx).genp;
      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);

      double angle = rFromNuMI.Angle(vec_p_mu);

      ret = TMath::Cos( angle );
    }

    return ret;
  });
  const TruthVar kTruth_MuonP([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_MuonIndex(nu);
    if(truth_idx>=0){
      const auto& mu_p = nu->prim.at(truth_idx).genp;
      ret = sqrt(mu_p.x*mu_p.x + mu_p.y*mu_p.y + mu_p.z*mu_p.z);
    }

    return ret;
  });
  const TruthVar kTruth_MuonPt([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_MuonIndex(nu);
    if(truth_idx>=0){
      const auto& p_mu = nu->prim.at(truth_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit(); // make unit vector

      double p_l = vec_p_mu.Dot(vec_p_nu);
      TVector3 vec_p_t_mu = vec_p_mu - p_l*vec_p_nu;

      ret = vec_p_t_mu.Mag();

    }

    return ret;
  });
  const TruthVar kTruth_MuonKE([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_MuonIndex(nu);
    if(truth_idx>=0){
      ret = nu->prim.at(truth_idx).genE - M_MUON;
    }

    return ret;
  });
  const TruthVar kTruth_MuonLength([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_MuonIndex(nu);
    if(truth_idx>=0){
      ret = nu->prim.at(truth_idx).length;
    }

    return ret;
  });
  const TruthVar kTruth_MuonContained([](const caf::SRTrueInteractionProxy *nu) -> int {
    // 1: contained, 0: not contained (-1: muon not found)
    int ret = -1;

    int truth_idx = kTruth_MuonIndex(nu);
    if(truth_idx>=0){
      if(nu->prim.at(truth_idx).contained) ret = 1;
      else ret = 0;
    }

    return ret;
  });

  // (Leading) Proton

  const TruthVar kTruth_ProtonIndex([](const caf::SRTrueInteractionProxy *nu) -> int {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < nu->prim.size(); ++i){
      // primary
      if( nu->prim.at(i).start_process!=0 ) continue;
      // proton
      if( abs(nu->prim.at(i).pdg)!=2212 ) continue;
      // non-nan genE
      if(isnan(nu->prim.at(i).genE)) continue;

      double this_E = nu->prim.at(i).genE;
      // if larger E, update
      if(this_E>max_E){
        max_E = this_E;
        truth_idx = i;
      }
    }
    return truth_idx;
  });
  const TruthVar kTruth_ProtonNuCosineTheta([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ProtonIndex(nu);
    if(truth_idx>=0){

      const auto& p_pro = nu->prim.at(truth_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      ret = vec_p_pro.Unit().Dot( vec_p_nu.Unit() );

    }

    return ret;
  });
  const TruthVar kTruth_ProtonCosThBeam([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ProtonIndex(nu);
    if(truth_idx>=0){

      TVector3 rFromNuMI(315.120380, 33.644912, 733.632532);

      const auto& p_pro = nu->prim.at(truth_idx).genp;
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);

      double angle = rFromNuMI.Angle(vec_p_pro);

      ret = TMath::Cos( angle );
    }

    return ret;
  });
  const TruthVar kTruth_ProtonP([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ProtonIndex(nu);
    if(truth_idx>=0){
      const auto& pro_p = nu->prim.at(truth_idx).genp;
      ret = sqrt(pro_p.x*pro_p.x + pro_p.y*pro_p.y + pro_p.z*pro_p.z);
    }

    return ret;
  });
  const TruthVar kTruth_ProtonPt([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ProtonIndex(nu);
    if(truth_idx>=0){
      const auto& p_pro = nu->prim.at(truth_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit();

      double p_l = vec_p_pro.Dot(vec_p_nu);
      TVector3 vec_p_t_pro = vec_p_pro - p_l*vec_p_nu;

      ret = vec_p_t_pro.Mag();

    }

    return ret;
  });
  const TruthVar kTruth_ProtonKE([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ProtonIndex(nu);
    if(truth_idx>=0){
      ret = nu->prim.at(truth_idx).genE - M_PROTON;
    }

    return ret;
  });
  const TruthVar kTruth_ProtonLength([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ProtonIndex(nu);
    if(truth_idx>=0){
      ret = nu->prim.at(truth_idx).length;
    }

    return ret;
  });

  // Muon+Proton
  const TruthVar kTruth_CosThMuonProton([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_mu_idx = kTruth_MuonIndex(nu);
    int truth_pro_idx = kTruth_ProtonIndex(nu);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = nu->prim.at(truth_mu_idx).genp;
      const auto& p_pro = nu->prim.at(truth_pro_idx).genp;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);

      ret = vec_p_mu.Unit().Dot( vec_p_pro.Unit() );

    }

    return ret;
  });

  // TKI
  // https://arxiv.org/abs/1910.08658
  double CalcTKI_deltaPT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu){

    TVector3 unit_vec_p_nu = vec_p_nu.Unit();

    // Get transverse momenta w.r.t. the neutrino direction
    TVector3 pt_mu = vec_p_mu - (vec_p_mu.Dot(unit_vec_p_nu))*unit_vec_p_nu ;
    TVector3 pt_pro = vec_p_pro - (vec_p_pro.Dot(unit_vec_p_nu))*unit_vec_p_nu;

    TVector3 vec_deltaPT = pt_mu+pt_pro;

    return vec_deltaPT.Mag();

  }
  double CalcTKI_deltaPTx(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu){

    TVector3 unit_vec_p_nu = vec_p_nu.Unit();

    // Get transverse momenta w.r.t. the neutrino direction
    TVector3 pt_mu = vec_p_mu - (vec_p_mu.Dot(unit_vec_p_nu))*unit_vec_p_nu ;
    TVector3 pt_pro = vec_p_pro - (vec_p_pro.Dot(unit_vec_p_nu))*unit_vec_p_nu;

    TVector3 vec_deltaPT = pt_mu+pt_pro;

    double deltaPT_x = ( unit_vec_p_nu.Cross(pt_mu.Unit()) ).Dot(vec_deltaPT);

    return deltaPT_x;

  }
  double CalcTKI_deltaPTy(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu){

    TVector3 unit_vec_p_nu = vec_p_nu.Unit();

    // Get transverse momenta w.r.t. the neutrino direction
    TVector3 pt_mu = vec_p_mu - (vec_p_mu.Dot(unit_vec_p_nu))*unit_vec_p_nu ;
    TVector3 pt_pro = vec_p_pro - (vec_p_pro.Dot(unit_vec_p_nu))*unit_vec_p_nu;

    TVector3 vec_deltaPT = pt_mu+pt_pro;

    double deltaPT_y = -1.*(pt_mu.Unit().Dot(vec_deltaPT));

    return deltaPT_y;

  }
  double CalcTKI_deltaalphaT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu){

    TVector3 unit_vec_p_nu = vec_p_nu.Unit();

    // Get transverse momenta w.r.t. the neutrino direction
    TVector3 pt_mu = vec_p_mu - (vec_p_mu.Dot(unit_vec_p_nu))*unit_vec_p_nu ;
    TVector3 pt_pro = vec_p_pro - (vec_p_pro.Dot(unit_vec_p_nu))*unit_vec_p_nu;

    TVector3 vec_deltaPT = pt_mu+pt_pro;

    double CosdeltaalphaT = -1. * pt_mu.Unit().Dot( vec_deltaPT.Unit() );
    double deltaalphaT = TMath::ACos( CosdeltaalphaT );
    return deltaalphaT;

  }
  double CalcTKI_deltaphiT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu){

    TVector3 unit_vec_p_nu = vec_p_nu.Unit();

    // Get transverse momenta w.r.t. the neutrino direction
    TVector3 pt_mu = vec_p_mu - (vec_p_mu.Dot(unit_vec_p_nu))*unit_vec_p_nu ;
    TVector3 pt_pro = vec_p_pro - (vec_p_pro.Dot(unit_vec_p_nu))*unit_vec_p_nu;

    TVector3 vec_deltaPT = pt_mu+pt_pro;

    double CosdeltaphiT = -1. * pt_mu.Unit().Dot( pt_pro.Unit() );
    double deltaphiT = TMath::ACos( CosdeltaphiT );
    return deltaphiT;
  }

  // TKI
  const TruthVar kTruth_deltaPT([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_mu_idx = kTruth_MuonIndex(nu);
    int truth_pro_idx = kTruth_ProtonIndex(nu);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = nu->prim.at(truth_mu_idx).genp;
      const auto& p_pro = nu->prim.at(truth_pro_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret =  CalcTKI_deltaPT(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  });
  const TruthVar kTruth_deltaPTx([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-99999); // TODO deltaPTx is a signed variable.. for now just making it very very largly negative

    int truth_mu_idx = kTruth_MuonIndex(nu);
    int truth_pro_idx = kTruth_ProtonIndex(nu);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = nu->prim.at(truth_mu_idx).genp;
      const auto& p_pro = nu->prim.at(truth_pro_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaPTx(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  });
  const TruthVar kTruth_deltaPTy([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-99999); // TODO deltaPTy is a signed variable.. for now just making it very very largly negative

    int truth_mu_idx = kTruth_MuonIndex(nu);
    int truth_pro_idx = kTruth_ProtonIndex(nu);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = nu->prim.at(truth_mu_idx).genp;
      const auto& p_pro = nu->prim.at(truth_pro_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaPTy(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  });
  const TruthVar kTruth_deltaalphaT([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f); // acos is in the interval of [0,pi]

    int truth_mu_idx = kTruth_MuonIndex(nu);
    int truth_pro_idx = kTruth_ProtonIndex(nu);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = nu->prim.at(truth_mu_idx).genp;
      const auto& p_pro = nu->prim.at(truth_pro_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaalphaT(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  });
  const TruthVar kTruth_deltaphiT([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f); // acos is in the interval of [0,pi]

    int truth_mu_idx = kTruth_MuonIndex(nu);
    int truth_pro_idx = kTruth_ProtonIndex(nu);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = nu->prim.at(truth_mu_idx).genp;
      const auto& p_pro = nu->prim.at(truth_pro_idx).genp;
      const auto& p_nu = nu->momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaphiT(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  });

  // (Leading) Charged pion

  const TruthVar kTruth_ChargedPionIndex([](const caf::SRTrueInteractionProxy *nu) -> int {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < nu->prim.size(); ++i){
      // primary
      if( nu->prim.at(i).start_process!=0 ) continue;
      // pi+-
      if( abs(nu->prim.at(i).pdg)!=211 ) continue;
      // non-nan genE
      if(isnan(nu->prim.at(i).genE)) continue;

      double this_E = nu->prim.at(i).genE;
      // if larger E, update
      if(this_E>max_E){
        max_E = this_E;
        truth_idx = i;
      }
    }
    return truth_idx;
  });
  const TruthVar kTruth_ChargedPionKE([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ChargedPionIndex(nu);
    if(truth_idx>=0){
      ret = nu->prim.at(truth_idx).genE - M_CHARGEDPION;
    }

    return ret;
  });

} // end namespace ana
