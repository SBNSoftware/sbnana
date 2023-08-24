#include "sbnana/SBNAna/Vars/PrimaryUtils.h"

using namespace std;
using namespace ana;

namespace ana{

namespace PrimaryUtil{

  // Interaction
  double NeutrinoE(const TrueInteraction& true_int){
    return true_int.E;
  }
  int NeutrinoPDG(const TrueInteraction& true_int){
    return true_int.pdg;
  }
  int NeutrinoMode(const TrueInteraction& true_int){
    return true_int.genie_mode;
  }
  int Target(const TrueInteraction& true_int){
    return true_int.targetPDG;
  }
  int Npip(const TrueInteraction& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 211 ) NPtl++;
    }
    return NPtl;
  }
  int Npim(const TrueInteraction& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == -211 ) NPtl++;
    }
    return NPtl;
  }
  int Npi0(const TrueInteraction& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( abs(prim.pdg) == 111 ) NPtl++;
    }
    return NPtl;
  }

  // Muon

  int MuonIndex(const TrueInteraction& true_int){
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

  double MuonNuCosineTheta(const TrueInteraction& true_int){
    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      return vec_p_mu.Unit().Dot( vec_p_nu.Unit() );

    }
    else{
      return -999.;
    }
  }
  double MuonP(const TrueInteraction& true_int){
    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){
      const auto& mu_p = true_int.prim.at(truth_idx).genp;
      return sqrt(mu_p.x*mu_p.x + mu_p.y*mu_p.y + mu_p.z*mu_p.z);
    }
    else{
      return -999.;
    }
  }
  double MuonPt(const TrueInteraction& true_int){
    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){
      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit(); // make unit vector

      double p_l = vec_p_mu.Dot(vec_p_nu);
      TVector3 vec_p_l_mu = vec_p_mu - p_l*vec_p_nu;

      return vec_p_l_mu.Mag();

    }
    else{
      return -999.;
    }
  }
  double MuonCosThBeam(const TrueInteraction& true_int){
    int truth_idx = MuonIndex(true_int);
    if(truth_idx>=0){

      TVector3 rFromNuMI(315.120380, 33.644912, 733.632532);

      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);

      double angle = rFromNuMI.Angle(vec_p_mu);

      return TMath::Cos( angle );
    }
    else{
      return -999.;
    }
  }

  // Proton

  int ProtonIndex(const TrueInteraction& true_int){
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

  double ProtonNuCosineTheta(const TrueInteraction& true_int){
    int truth_idx = ProtonIndex(true_int);
    if(truth_idx>=0){

      const auto& p_pro = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      return vec_p_pro.Unit().Dot( vec_p_nu.Unit() );

    }
    else{
      return -999.;
    }
  }
  double ProtonP(const TrueInteraction& true_int){
    int truth_idx = ProtonIndex(true_int);
    if(truth_idx>=0){
      const auto& pro_p = true_int.prim.at(truth_idx).genp;
      return sqrt(pro_p.x*pro_p.x + pro_p.y*pro_p.y + pro_p.z*pro_p.z);
    }
    else{
      return -999.;
    }
  }
  double ProtonPt(const TrueInteraction& true_int){
    int truth_idx = ProtonIndex(true_int);
    if(truth_idx>=0){
      const auto& p_pro = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit();

      double p_l = vec_p_pro.Dot(vec_p_nu);
      TVector3 vec_p_l_pro = vec_p_pro - p_l*vec_p_nu;

      return vec_p_l_pro.Mag();

    }
    else{
      return -999.;
    }
  }

  // Muon+Proton
  double CosThMuonProton(const TrueInteraction& true_int){
    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);

      return vec_p_mu.Unit().Dot( vec_p_pro.Unit() );

    }
    else{
      return -999.;
    }
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

  double deltaPT(const TrueInteraction& true_int){
    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      return CalcTKI_deltaPT(vec_p_mu, vec_p_pro, vec_p_n);

    }
    else{
      return -999.;
    }
  }
  double deltaPTx(const TrueInteraction& true_int){
    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      return CalcTKI_deltaPTx(vec_p_mu, vec_p_pro, vec_p_n);

    }
    else{
      return -999.;
    }
  }
  double deltaPTy(const TrueInteraction& true_int){
    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      return CalcTKI_deltaPTy(vec_p_mu, vec_p_pro, vec_p_n);

    }
    else{
      return -999.;
    }
  }
  double deltaalphaT(const TrueInteraction& true_int){
    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      return CalcTKI_deltaalphaT(vec_p_mu, vec_p_pro, vec_p_n);

    }
    else{
      return -999.;
    }
  }
  double deltaphiT(const TrueInteraction& true_int){
    int truth_mu_idx = MuonIndex(true_int);
    int truth_pro_idx = ProtonIndex(true_int);
    if(truth_mu_idx>=0 && truth_pro_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_mu_idx).genp;
      const auto& p_pro = true_int.prim.at(truth_pro_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_n(p_nu.x, p_nu.y, p_nu.z);

      return CalcTKI_deltaphiT(vec_p_mu, vec_p_pro, vec_p_n);

    }
    else{
      return -999.;
    }
  }

} // END namespace PrimaryUtil

} // ENd namespace ana
