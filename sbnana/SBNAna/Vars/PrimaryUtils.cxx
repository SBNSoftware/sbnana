#include "sbnana/SBNAna/Vars/PrimaryUtils.h"

namespace ana{

namespace PrimaryUtil{

  // Interaction
  double NeutrinoE(const caf::SRTrueInteractionProxy& true_int){
    return true_int.E;
  }
  int NeutrinoPDG(const caf::SRTrueInteractionProxy& true_int){
    return true_int.pdg;
  }
  int NeutrinoMode(const caf::SRTrueInteractionProxy& true_int){
    return true_int.genie_mode;
  }
  int Target(const caf::SRTrueInteractionProxy& true_int){
    return true_int.targetPDG;
  }
  int Npip(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 211 ) NPtl++;
    }
    return NPtl;
  }
  int Npim(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == -211 ) NPtl++;
    }
    return NPtl;
  }
  int Npi0(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( abs(prim.pdg) == 111 ) NPtl++;
    }
    return NPtl;
  }

  // Muon

  int MuonIndex(const caf::SRTrueInteractionProxy& true_int){
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < true_int.prim.size(); ++i){
      // primary
      if( true_int.prim.at(i).start_process!=0 ) continue;
      // muon
      if( abs(true_int.prim.at(i).pdg)!=13 ) continue;
      // non-nan genE
      if(isnan(true_int.prim.at(i).genE)) continue;

      double this_E = true_int.prim.at(i).genE;
      // if larger E, update
      if(this_E>max_E){
        max_E = this_E;
        truth_idx = i;
      }
    }
    return truth_idx;
  }

  double MuonNuCosineTheta(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      ret = vec_p_mu.Unit().Dot( vec_p_nu.Unit() );

    }

    return ret;
  }
  double MuonP(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){
      const auto& mu_p = true_int.prim.at(truth_idx).genp;
      ret = sqrt(mu_p.x*mu_p.x + mu_p.y*mu_p.y + mu_p.z*mu_p.z);
    }

    return ret;

  }
  double MuonPt(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){
      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit(); // make unit vector

      double p_l = vec_p_mu.Dot(vec_p_nu);
      TVector3 vec_p_l_mu = vec_p_mu - p_l*vec_p_nu;

      ret = vec_p_l_mu.Mag();

    }

    return ret;

  }
  double MuonCosThBeam(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){

      TVector3 rFromNuMI(315.120380, 33.644912, 733.632532);

      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);

      double angle = rFromNuMI.Angle(vec_p_mu);

      ret = TMath::Cos( angle );
    }

    return ret;
  }

  // Proton

  int ProtonIndex(const caf::SRTrueInteractionProxy& true_int){
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < true_int.prim.size(); ++i){
      // primary
      if( true_int.prim.at(i).start_process!=0 ) continue;
      // proton
      if( abs(true_int.prim.at(i).pdg)!=2212 ) continue;
      // non-nan genE
      if(isnan(true_int.prim.at(i).genE)) continue;

      double this_E = true_int.prim.at(i).genE;
      // if larger E, update
      if(this_E>max_E){
        max_E = this_E;
        truth_idx = i;
      }
    }
    return truth_idx;
  }

  double ProtonNuCosineTheta(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = ProtonIndex(true_int);
    if(truth_idx>=0){

      const auto& p_pro = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      ret = vec_p_pro.Unit().Dot( vec_p_nu.Unit() );

    }

    return ret;
  }
  double ProtonP(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = ProtonIndex(true_int);
    if(truth_idx>=0){
      const auto& pro_p = true_int.prim.at(truth_idx).genp;
      ret = sqrt(pro_p.x*pro_p.x + pro_p.y*pro_p.y + pro_p.z*pro_p.z);
    }

    return ret;
  }
  double ProtonPt(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = ProtonIndex(true_int);
    if(truth_idx>=0){
      const auto& p_pro = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit();

      double p_l = vec_p_pro.Dot(vec_p_nu);
      TVector3 vec_p_l_pro = vec_p_pro - p_l*vec_p_nu;

      ret = vec_p_l_pro.Mag();

    }

    return ret;
  }

  // Muon+Proton
  double CosThMuonProton(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);

      ret = vec_p_mu.Unit().Dot( vec_p_pro.Unit() );

    }

    return ret;
  }

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

  double deltaPT(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret =  CalcTKI_deltaPT(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  }
  double deltaPTx(const caf::SRTrueInteractionProxy& true_int){

    double ret(-99999); // TODO deltaPTx is a signed variable.. for now just making it very very largly negative

    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaPTx(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  }
  double deltaPTy(const caf::SRTrueInteractionProxy& true_int){

    double ret(-99999); // TODO deltaPTy is a signed variable.. for now just making it very very largly negative

    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaPTy(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  }
  double deltaalphaT(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f); // acos is in the interval of [0,pi]

    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaalphaT(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  }
  double deltaphiT(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f); // acos is in the interval of [0,pi]

    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      ret = CalcTKI_deltaphiT(vec_p_mu, vec_p_pro, vec_p_n);

    }

    return ret;
  }

} // END namespace PrimaryUtil

} // ENd namespace ana
