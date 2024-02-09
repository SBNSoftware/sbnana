#include "sbnana/SBNAna/Cuts/NuMIXSecSysts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecTruthVars.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <iostream>

namespace ana {

  // NuMIXSecPiSyst

  NuMIXSecPiSyst::NuMIXSecPiSyst(const std::string& name, const std::string& latexName):
    ISyst(name, latexName)
  {

    kPdgNeutron = 2112;
    kPdgProton = 2212;
    kPdgClusterNN = 2000000200;
    kPdgClusterNP = 2000000201;
    kPdgClusterPP = 2000000202;

    kPdgPiP =   211; // pi+
    kPdgPiM =  -211; // pi-
    kPdgPi0 =   111; // pi0
    kPdgEta =   221; // eta
    kPdgKP =   321; // K+
    kPdgKM =  -321; // K-
    kPdgK0 =   311; // K0
    kPdgAntiK0 =  -311; // \bar{K0}
    kPdgLambda =  3122; // Lambda
    kPdgAntiLambda = -3122; // \bar{Lambda}
    kPdgGamma =    22; // photon

  }

  double NuMIXSecPiSyst::GetSPPQ2Reweight(double Q2_GeV2, double parameter_value) const {

    double X = Q2_GeV2;
    if(X>=3.000000) X = 3.000000;

    double this_rw = 1.;
    if( X < 0.025000) this_rw = 1.253255;
    else if( X >= 0.025000 && X < 0.050000) this_rw = 1.589738;
    else if( X >= 0.050000 && X < 0.100000) this_rw = 1.733869;
    else if( X >= 0.100000 && X < 0.200000) this_rw = 1.651728;
    else if( X >= 0.200000 && X < 0.300000) this_rw = 1.659705;
    else if( X >= 0.300000 && X < 0.400000) this_rw = 1.584229;
    else if( X >= 0.400000 && X < 0.500000) this_rw = 1.703793;
    else if( X >= 0.500000 && X < 0.700000) this_rw = 1.475510;
    else if( X >= 0.700000 && X < 1.000000) this_rw = 1.456727;
    else if( X >= 1.000000 && X < 1.300000) this_rw = 1.252215;
    else if( X >= 1.300000 && X < 2.000000) this_rw = 1.048199;
    else if( X >= 2.000000 && X < 3.000000) this_rw = 1.650489;
    else{
      this_rw = 1.;
    }

    return (1-parameter_value) * 1. + parameter_value * this_rw;

  }
  double NuMIXSecPiSyst::GetSPPTpiReweight(double Tpi_GeV, double parameter_value) const {

    static double const P0 = 1.347965;
    //static double const P0Err = 0.097783;
    static double const P1 = -2.817372;
    //static double const P1Err = 0.320814;

    double X = Tpi_GeV;
    if(X>0.350) X = 0.350;

    double this_rw = P0 + P1 * X;
    if(this_rw<0) this_rw = 1.;

    //std::cout << "[GetSPPTpiReweight] Tpi_GeV = " << Tpi_GeV << ", RW = " << this_rw << std::endl;

    // dial = 0 : 1
    // dial = 1 : this_rw

    return (1-parameter_value) * 1. + parameter_value * this_rw;

  }

  void NuMIXSecPiSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    this->Shift(sigma, &sr->truth, weight);
  }

  void NuMIXSecPiSyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

    // Check res

    int NeutCode = NeutReactionCode(sr);
    if( abs(NeutCode)!=11 || abs(NeutCode)!= 12 || abs(NeutCode)!= 13) return;

    int TargetPDG = kTruth_Target(sr);
    int TargetA = ((TargetPDG % 10000)) / 10;
    if(TargetA==1) return;

    double Q2 = kTruth_Q2(sr);
    double Tpi = kTruth_ChargedPionKE(sr);

    double Q2RW = GetSPPQ2Reweight(Q2, sigma);
    double TpiRW = GetSPPTpiReweight(Tpi, sigma);

    weight *= Q2RW*TpiRW;

  }

  int NuMIXSecPiSyst::NeutReactionCode(caf::SRTrueInteractionProxy *sr) const {

    int evtype = 0;

    int GENIEMode = kTruth_NeutrinoMode(sr);
    int GENIEIntType = sr->genie_inttype;
    int HitNucPDG = sr->hitnuc;
    int NuPDG = sr->pdg;

    bool is_cc    = sr->iscc;
    bool is_nc    = sr->isnc;
    bool is_charm = sr->ischarm;

    bool is_qel   = GENIEMode==caf::kQE;
    bool is_dis   = GENIEMode==caf::kDIS;
    bool is_res   = GENIEMode==caf::kRes;
    bool is_coh_pr = GENIEMode==caf::kCoh;
    bool is_ve    = GENIEIntType==caf::kNuElectronElastic;
    bool is_mec   = GENIEMode==caf::kMEC;
    bool is_imd   = GENIEIntType==caf::kInverseMuDecay;

    //bool is_ask   = proc.IsSingleKaon(); 
    bool is_ask = false;

    bool is_diff  = GENIEMode==caf::kDiffractive;

    bool is_p     = HitNucIsSet(HitNucPDG) ? HitNucPDG==kPdgProton  : false;
    bool is_n     = HitNucIsSet(HitNucPDG) ? HitNucPDG==kPdgNeutron : false;

    bool is_nu    = NuPDG>0;
    bool is_nubar = NuPDG<0;

    bool W_gt_2 = ( isnan(sr->w) || isinf(sr->w) ) ? false : sr->w>2.0;

    // (quasi-)elastic, nc+cc, nu+nubar
    //
    if      (is_qel && !is_charm && is_cc && is_nu           ) evtype =   1;
    else if (is_qel && !is_charm && is_nc && is_nu && is_p   ) evtype =  51;
    else if (is_qel && !is_charm && is_nc && is_nu && is_n   ) evtype =  52;
    else if (is_qel && !is_charm && is_cc && is_nubar        ) evtype =  -1;
    else if (is_qel && !is_charm && is_nc && is_nubar && is_p) evtype = -51;
    else if (is_qel && !is_charm && is_nc && is_nubar && is_n) evtype = -52;

    // MEC - only CC implemented in NEUT
    else if      (is_mec && !is_charm && is_cc && is_nu           ) evtype =   2;
    else if      (is_mec && !is_charm && is_cc && is_nubar        ) evtype =  -2;

    // quasi-elastic charm production
    // Part of the DIS W>2GeV mode in NEUT - CB
    else if (is_qel && is_charm && is_cc && is_nu    ) evtype =   26;
    else if (is_qel && is_charm && is_cc && is_nubar ) evtype =  -26;

    // inverse mu- (tau-) decay and ve- elastic
    //Those modes don't actually exist in NEUT, 9 and 59 used as place holders
    else if ( is_imd ) evtype =  9;
    else if ( is_ve  ) evtype = 59;

    // coherent pi, nc+cc, nu+nubar
    //
    else if (is_coh_pr && is_cc && is_nu   ) evtype =  16;
    else if (is_coh_pr && is_cc && is_nubar) evtype = -16;
    else if (is_coh_pr && is_nc && is_nu   ) evtype =  36;
    else if (is_coh_pr && is_nc && is_nubar) evtype = -36;

    // dis, W>2, nc+cc, nu+nubar
    // (charm DIS not simulated by NEUT, will bundle GENIE charm DIS into this category)
    //
    else if (is_dis && W_gt_2 && is_cc && is_nu   ) evtype =  26;
    else if (is_dis && W_gt_2 && is_nc && is_nu   ) evtype =  46;
    else if (is_dis && W_gt_2 && is_cc && is_nubar) evtype = -26;
    else if (is_dis && W_gt_2 && is_nc && is_nubar) evtype = -46;

    // resonance or dis with W < 2 GeV or single kaon
    //
    else if ( is_res || (is_dis && !W_gt_2) || is_ask ) {

      int nn=0, np=0, npi0=0, npip=0, npim=0, nKp=0, nKm=0, nK0=0, neta=0, nlambda=0, ngamma=0;

      for(const auto& prim: sr->prim){

        // only primary
        if ( prim.start_process != 0 ) continue;

        int ghep_pdgc = prim.pdg;

        // We can't replicate "bool count_it" part from caf primaries..

        if(ghep_pdgc == kPdgProton )    np++;            // p
        if(ghep_pdgc == kPdgNeutron)    nn++;            // n
        if(ghep_pdgc == kPdgPiP)        npip++;          // pi+
        if(ghep_pdgc == kPdgPiM)        npim++;          // pi-
        if(ghep_pdgc == kPdgPi0)        npi0++;          // pi0
        if(ghep_pdgc == kPdgEta)        neta++;          // eta0
        if(ghep_pdgc == kPdgKP)         nKp++;           // K+
        if(ghep_pdgc == kPdgKM)         nKm++;           // K-
        if(ghep_pdgc == kPdgK0)         nK0++;           // K0
        if(ghep_pdgc == kPdgAntiK0)     nK0++;           // K0
        if(ghep_pdgc == kPdgLambda)     nlambda++;       // Lamda
        if(ghep_pdgc == kPdgAntiLambda) nlambda++;       // Lamda
        if(ghep_pdgc == kPdgGamma)      ngamma++;        // photon

      }
    
      int nnuc = np + nn;
      int npi  = npi0 + npip + npim;
      int nK   = nK0 + nKp + nKm;
      int neKL = neta + nK + nlambda;

      bool is_radiative_dec = (nnuc==1) && (npi==0) && (ngamma==1);

      //
      // single gamma from resonances
      //

      if      (is_res && is_nu    && is_cc && is_n && is_radiative_dec) evtype =  17;
      else if (is_res && is_nu    && is_nc && is_n && is_radiative_dec) evtype =  38;
      else if (is_res && is_nu    && is_nc && is_p && is_radiative_dec) evtype =  39;
      
      else if (is_res && is_nubar && is_cc && is_p && is_radiative_dec) evtype = -17;
      else if (is_res && is_nubar && is_nc && is_n && is_radiative_dec) evtype = -38;
      else if (is_res && is_nubar && is_nc && is_p && is_radiative_dec) evtype = -39;
      
      //
      // single pi (res + non-res bkg)
      //
      
      // nu CC
      else if (is_nu    && is_cc && is_p && np==1 && nn==0 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  11;
      else if (is_nu    && is_cc && is_n && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  12;
      else if (is_nu    && is_cc && is_n && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  13;
      
      // nu NC
      else if (is_nu    && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  31;
      else if (is_nu    && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  32;
      else if (is_nu    && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype =  33;
      else if (is_nu    && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  34;
      
      //nubar CC
      else if (is_nubar && is_cc && is_n && np==0 && nn==1 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -11;
      else if (is_nubar && is_cc && is_p && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -12;
      else if (is_nubar && is_cc && is_p && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -13;
      
      //nubar NC
      else if (is_nubar && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -31;
      else if (is_nubar && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -32;
      else if (is_nubar && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -33;
      else if (is_nubar && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype = -34;
      
      //
      // single eta from res
      //
      
      else if (is_res &&  is_nu    && is_cc && is_n && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  22;
      else if (is_res &&  is_nu    && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  42;
      else if (is_res &&  is_nu    && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  43;
      
      else if (is_res &&  is_nubar && is_cc && is_p && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -22;
      else if (is_res &&  is_nubar && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -42;
      else if (is_res &&  is_nubar && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -43;
      
      //
      // single K from res (dS=0)
      //
      
      else if (is_res &&  is_nu    && is_cc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  23;
      else if (is_res &&  is_nu    && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  44;
      else if (is_res &&  is_nu    && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  45;
      
      else if (is_res &&  is_nubar && is_cc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -23;
      else if (is_res &&  is_nubar && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -44;
      else if (is_res &&  is_nubar && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -45;
      
      //
      // single K from AtharSingleKaon (dS=1)
      //
      //Those modes are assigned but not used (xsec=0) in NEUT
      else if (is_ask &&  is_nu && is_cc && is_n && nn==1 && np==0 && nKp==1 && neKL==1) evtype =  18;
      else if (is_ask &&  is_nu && is_cc && is_n && nn==0 && np==1 && nK0==1 && neKL==1) evtype =  19;
      else if (is_ask &&  is_nu && is_cc && is_p && nn==0 && np==1 && nKp==1 && neKL==1) evtype =  20;
      
      
      // antineutrino modes not yet implemented
      //else if (is_ask &&  is_nubar && is_cc && is_n && nn==1 && np==0 && nKp==1 && neKL==1) evtype = -18;
      //else if (is_ask &&  is_nubar && is_cc && is_n && nn==0 && np==1 && nK0==1 && neKL==1) evtype = -19;
      //else if (is_ask &&  is_nubar && is_cc && is_p && nn==0 && np==1 && nKp==1 && neKL==1) evtype = -20;
      
      //
      // multi-pi (res or dis (W<2GeV)
      //
      
      else if (is_nu    && is_cc && npi>1) evtype =  21;
      else if (is_nu    && is_nc && npi>1) evtype =  41;
      else if (is_nubar && is_cc && npi>1) evtype = -21;
      else if (is_nubar && is_nc && npi>1) evtype = -41;
      
      //
      // rare final state for RES or low-W (<2GeV) DIS events
      // (eg K0\bar{K0} final states, N0(1720) -> Sigma- K+ res decays, etc)
      // bundled-in with multi-pi
      //
      else {
        if (is_nu    && is_cc) evtype =  21;
        else if (is_nu    && is_nc) evtype =  41;
        else if (is_nubar && is_cc) evtype = -21;
        else if (is_nubar && is_nc) evtype = -41;
      }

    } // END red or dis
      
    // Weak diffractive processes
    else if ( is_diff && is_cc ) {
      if ( is_nu ) evtype = 15;
      else if ( is_nubar ) evtype = -15;
    }
      else if ( is_diff && is_nc ) {
      if ( is_nu ) evtype = 35;
      else if ( is_nubar ) evtype = -35;
    }
    
    return evtype;
      
  }

  bool NuMIXSecPiSyst::HitNucIsSet(int pdgc) const {
    bool ok =
    IsNucleon(pdgc)          ||
    Is2NucleonCluster(pdgc);
    return ok;
  }
  bool NuMIXSecPiSyst::IsNucleon(int pdgc) const {
    return (pdgc == kPdgNeutron || pdgc == kPdgProton);
  }
  bool NuMIXSecPiSyst::Is2NucleonCluster(int pdgc) const {
    return (
      pdgc == kPdgClusterNN   ||
      pdgc == kPdgClusterNP   ||
      pdgc == kPdgClusterPP
    );
  }

} // end namespace ana
