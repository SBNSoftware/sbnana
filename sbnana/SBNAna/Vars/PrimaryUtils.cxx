#include "sbnana/SBNAna/Vars/PrimaryUtils.h"

namespace ana{

namespace PrimaryUtil{

  // Interaction
  double NeutrinoE_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.E;
  }
  int NeutrinoPDG_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.pdg;
  }
  int NeutrinoMode_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.genie_mode;
  }
  int IsCC_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.iscc;
  }
  int Target_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.targetPDG;
  }

  int NMuon_True(const caf::SRTrueInteractionProxy& true_int){ //slc->truth.prim
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 13 ) NPtl++;
    }
    return NPtl;
  }

  int NProton_True(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 2212 ) NPtl++;
    }
    return NPtl;
  }
  int NNeutron_True(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 2112 ) NPtl++;
    }
    return NPtl;
  }
  int Npip_True(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == 211 ) NPtl++;
    }
    return NPtl;
  }
  int Npim_True(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( prim.pdg == -211 ) NPtl++;
    }
    return NPtl;
  }
  int Npi0_True(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;
      if ( abs(prim.pdg) == 111 ) NPtl++;
    }
    return NPtl;
  }
  int Npi0_True_Any(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( abs(prim.pdg) == 111 ) NPtl++;
    }
    return NPtl;
  }
  int NPhoton_True(const caf::SRTrueInteractionProxy& true_int){
    int NPtl = 0;
    for ( auto const& prim : true_int.prim ) {
      if ( abs(prim.pdg) == 22 && prim.start_process == 3 ) NPtl++;
    }
    return NPtl;
  }
  double Q2_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.Q2;
  }
  double q0_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.q0_lab;
  }
  double q3_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.modq_lab;
  }
  double w_True(const caf::SRTrueInteractionProxy& true_int){
    return true_int.w;
  }

  // Muon

  int MuonIndex_True(const caf::SRTrueInteractionProxy& true_int){
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

  // True Neutral Pion

  int Pi0_G4ID(const caf::SRTrueInteractionProxy& true_int){
    double max_E(-999);
    int g4id(-1);
    for(std::size_t i(0); i < true_int.prim.size(); ++i){
      // primary
      if( true_int.prim.at(i).start_process!=0 ) continue;
      // it decays
      if( true_int.prim.at(i).end_process!=3 ) continue;
      // pi0
      if( abs(true_int.prim.at(i).pdg) != 111 ) continue;
      // non-nan genE
      if(isnan(true_int.prim.at(i).genE)) continue;
      // daughters
      if(true_int.prim.at(i).daughters.size() != 2) continue;

      double this_E = true_int.prim.at(i).genE;
      // if larger E, update
      if(this_E>max_E){
        max_E = this_E;
        g4id = true_int.prim.at(i).G4ID;
      }
    }
    return g4id;
  }

  int Pi0LeadingPhotonG4ID(const caf::SRTrueInteractionProxy& true_int){
    int parentg4id = Pi0_G4ID(true_int);
    if(parentg4id<0) return -1;
    int leadingg4id = -1;
    double leadingE(-999);
    for(std::size_t i(0); i < true_int.prim.size(); ++i){
      // Is Photon
      if( abs(true_int.prim.at(i).pdg) != 22 ) continue;
      // non-nan genE
      if(isnan(true_int.prim.at(i).genE)) continue;
      // parent G4ID
      if( (int) true_int.prim.at(i).parent != parentg4id) continue;
      double this_E = true_int.prim.at(i).genE;
      // if larger E, update
      if(this_E>leadingE){
        leadingE = this_E;
        leadingg4id = true_int.prim.at(i).G4ID;
      }
    }
    return leadingg4id;
  }

  int Pi0SubLeadingPhotonG4ID(const caf::SRTrueInteractionProxy& true_int){
    int parentg4id = Pi0_G4ID(true_int);
    if(parentg4id<0) return -1;
    int subleadingg4id = -1;
    double subleadingE(-999);
    for(std::size_t i(0); i < true_int.prim.size(); ++i){
      // Is Photon
      if( abs(true_int.prim.at(i).pdg) != 22 ) continue;
      // non-nan genE
      if(isnan(true_int.prim.at(i).genE)) continue;
      // parent G4ID
      if( (int) true_int.prim.at(i).parent != parentg4id) continue;
      // leading photon
      if(true_int.prim.at(i).G4ID == Pi0LeadingPhotonG4ID(true_int)) continue;
      double this_E = true_int.prim.at(i).genE;
      // if larger E, update
      if(this_E>subleadingE){
        subleadingE = this_E;
        subleadingg4id = true_int.prim.at(i).G4ID;
      }
    }
    return subleadingg4id;
  }

  double MuonNuCosineTheta_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex_True(true_int);
    if(truth_idx>=0){

      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      ret = vec_p_mu.Unit().Dot( vec_p_nu.Unit() );

    }

    return ret;
  }
  double MuonCosThBeam_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex_True(true_int);
    if(truth_idx>=0){

      TVector3 rFromNuMI(315.120380, 33.644912, 733.632532);

      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);

      double angle = rFromNuMI.Angle(vec_p_mu);

      ret = TMath::Cos( angle );
    }

    return ret;
  }
  double MuonP_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex_True(true_int);
    if(truth_idx>=0){
      const auto& mu_p = true_int.prim.at(truth_idx).genp;
      ret = sqrt(mu_p.x*mu_p.x + mu_p.y*mu_p.y + mu_p.z*mu_p.z);
    }

    return ret;

  }
  double MuonPt_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = MuonIndex_True(true_int);
    if(truth_idx>=0){
      const auto& p_mu = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_mu(p_mu.x, p_mu.y, p_mu.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit(); // make unit vector

      double p_l = vec_p_mu.Dot(vec_p_nu);
      TVector3 vec_p_t_mu = vec_p_mu - p_l*vec_p_nu;

      ret = vec_p_t_mu.Mag();

    }

    return ret;
  }
  double MuonKE_True(const caf::SRTrueInteractionProxy& true_int){
    double ret(-5.f);
    
    int truth_idx = MuonIndex_True(true_int);
    if(truth_idx>=0){
      ret = true_int.prim.at(truth_idx).genE - M_MUON;
    }
    
    return ret;
  }
  double MuonLength_True(const caf::SRTrueInteractionProxy& true_int){
    double ret(-5.f);

    int truth_idx = MuonIndex_True(true_int);
    if(truth_idx>=0){
      ret = true_int.prim.at(truth_idx).length;
    }

    return ret;
  }

  // Proton

  int ProtonIndex_True(const caf::SRTrueInteractionProxy& true_int){
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

  double ProtonNuCosineTheta_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = ProtonIndex_True(true_int);
    if(truth_idx>=0){

      const auto& p_pro = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);

      ret = vec_p_pro.Unit().Dot( vec_p_nu.Unit() );

    }

    return ret;
  }
  double ProtonCosThBeam_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = ProtonIndex_True(true_int);
    if(truth_idx>=0){

      TVector3 rFromNuMI(315.120380, 33.644912, 733.632532);

      const auto& p_pro = true_int.prim.at(truth_idx).genp;
      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);

      double angle = rFromNuMI.Angle(vec_p_pro);

      ret = TMath::Cos( angle );
    }

    return ret;
  }
  double ProtonP_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = ProtonIndex_True(true_int);
    if(truth_idx>=0){
      const auto& pro_p = true_int.prim.at(truth_idx).genp;
      ret = sqrt(pro_p.x*pro_p.x + pro_p.y*pro_p.y + pro_p.z*pro_p.z);
    }

    return ret;
  }
  double ProtonPt_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_idx = ProtonIndex_True(true_int);
    if(truth_idx>=0){
      const auto& p_pro = true_int.prim.at(truth_idx).genp;
      const auto& p_nu = true_int.momentum;

      TVector3 vec_p_pro(p_pro.x, p_pro.y, p_pro.z);
      TVector3 vec_p_nu(p_nu.x, p_nu.y, p_nu.z);
      vec_p_nu = vec_p_nu.Unit();

      double p_l = vec_p_pro.Dot(vec_p_nu);
      TVector3 vec_p_t_pro = vec_p_pro - p_l*vec_p_nu;

      ret = vec_p_t_pro.Mag();

    }

    return ret;
  }
  double ProtonKE_True(const caf::SRTrueInteractionProxy& true_int){
    double ret(-5.f);

    int truth_idx = ProtonIndex_True(true_int);
    if(truth_idx>=0){
      ret = true_int.prim.at(truth_idx).genE - M_PROTON;
    }

    return ret;
  }
  double ProtonLength_True(const caf::SRTrueInteractionProxy& true_int){
    double ret(-5.f);

    int truth_idx = ProtonIndex_True(true_int);
    if(truth_idx>=0){
      ret = true_int.prim.at(truth_idx).length;
    }

    return ret;
  }

  // Muon+Proton
  double CosThMuonProton_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_mu_idx = MuonIndex_True(true_int);
    int truth_pro_idx = ProtonIndex_True(true_int);
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

  double deltaPT_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f);

    int truth_mu_idx = MuonIndex_True(true_int);
    int truth_pro_idx = ProtonIndex_True(true_int);
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
  double deltaPTx_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-99999); // TODO deltaPTx is a signed variable.. for now just making it very very largly negative

    int truth_mu_idx = MuonIndex_True(true_int);
    int truth_pro_idx = ProtonIndex_True(true_int);
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
  double deltaPTy_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-99999); // TODO deltaPTy is a signed variable.. for now just making it very very largly negative

    int truth_mu_idx = MuonIndex_True(true_int);
    int truth_pro_idx = ProtonIndex_True(true_int);
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
  double deltaalphaT_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f); // acos is in the interval of [0,pi]

    int truth_mu_idx = MuonIndex_True(true_int);
    int truth_pro_idx = ProtonIndex_True(true_int);
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
  double deltaphiT_True(const caf::SRTrueInteractionProxy& true_int){

    double ret(-5.f); // acos is in the interval of [0,pi]

    int truth_mu_idx = MuonIndex_True(true_int);
    int truth_pro_idx = ProtonIndex_True(true_int);
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

  // Pi0
  // Truth angle between photons from neutral pion decay
  double CosThPhotonPhoton_True(const caf::SRTrueInteractionProxy& true_int){
    if(Pi0_G4ID(true_int)<0) return -999.;
    if(Pi0LeadingPhotonG4ID(true_int)<0) return -999.;
    if(Pi0SubLeadingPhotonG4ID(true_int)<0) return -999.;
    int leadingg4id = Pi0LeadingPhotonG4ID(true_int);
    int subleadingg4id = Pi0SubLeadingPhotonG4ID(true_int);
    double ret(-999.f);
    TVector3 leadingp, subleadingp;
    for(std::size_t i(0); i < true_int.prim.size(); ++i){
      if( (int) true_int.prim.at(i).G4ID == leadingg4id ){
        leadingp.SetXYZ(true_int.prim.at(i).end.x - true_int.prim.at(i).start.x, 
                        true_int.prim.at(i).end.y - true_int.prim.at(i).start.y,
                        true_int.prim.at(i).end.z - true_int.prim.at(i).start.z);
      }
      if( (int) true_int.prim.at(i).G4ID == subleadingg4id ){
        subleadingp.SetXYZ(true_int.prim.at(i).end.x - true_int.prim.at(i).start.x, 
                           true_int.prim.at(i).end.y - true_int.prim.at(i).start.y,
                           true_int.prim.at(i).end.z - true_int.prim.at(i).start.z);
      } 
    }
    if(leadingp.Mag() > 0 && subleadingp.Mag() > 0){
      ret = leadingp.Unit().Dot(subleadingp.Unit());
    }
    return ret;
  }

} // END namespace PrimaryUtil

} // ENd namespace ana
