#include "sbnana/SBNAna/Cuts/NuMIXSecSysts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecTruthVars.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <iostream>

namespace ana {

  // NuMIXSecPiSyst

  NuMIXSecPiSyst::NuMIXSecPiSyst(const std::string& name, const std::string& latexName):
    ISyst(name, latexName)
  {

  }

  double NuMIXSecPiSyst::GetSPPQ2Reweight(double Q2_GeV2) const {

    double X = Q2_GeV2;
    if(Q2_GeV2>=3.0) X = 2.5;

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
      this_rw = 1.650489;
    }

    return this_rw;

  }
  double NuMIXSecPiSyst::GetSPPTpiReweight(double Tpi_GeV) const {

    static double const P0 = 1.347965;
    //static double const P0Err = 0.097783;
    static double const P1 = -2.817372;
    //static double const P1Err = 0.320814;

    double X = Tpi_GeV;
    if(X>=0.350) X = 0.350;

    double this_rw = P0 + P1 * X;
    if(this_rw<0) this_rw = 1.;

    return this_rw;

  }

  void NuMIXSecPiSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    this->Shift(sigma, &sr->truth, weight);
  }

  void NuMIXSecPiSyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

    // Check process

    int genie_n_photons = 0;
    int genie_n_mesons = 0;
    int nPip = 0;
    for(const auto& prim: sr->prim){

      // only primary
      if ( prim.start_process != 0 ) continue;

      int pdg = prim.pdg;
      int energy = prim.genE * 1000.; // GeV->MeV

      if (pdg == 22 && energy > 10.0) {
        genie_n_photons++;
      }
      if (abs(pdg) == 211 || //pi+-
               pdg == 111 ||  // pi0
               abs(pdg) == 321 || // K-
               abs(pdg) == 323 || // K*+-
               pdg == 130 || // KL0
               pdg == 310 || // KS0
               pdg == 311 || // K0
               pdg == 313 || // K*0
               abs(pdg) == 221 || // eta
               abs(pdg) == 331 // eta' (958)
               ) {
        genie_n_mesons++;
      }
      if(pdg==211){
        nPip++;
      }

    }

    if (nPip != 1 || genie_n_mesons!= 1)
      return;
    if (genie_n_photons != 0 )
      return;


    int TargetPDG = kTruth_Target(sr);
    int TargetA = ((TargetPDG % 10000)) / 10;
    if(TargetA==1) return;

    double Q2 = kTruth_Q2(sr);
    double Tpi = kTruth_ChargedPionKE(sr);

    double Q2RW = GetSPPQ2Reweight(Q2);
    double TpiRW = GetSPPTpiReweight(Tpi);
    double FullRW = Q2RW*TpiRW;

    double this_sigma = sigma;
    if(sigma<0) this_sigma = 0.;
    if(sigma>1) this_sigma = 1.;

    double this_rw = (1.-this_sigma) * 1. + this_sigma * FullRW;

    weight *= this_rw;

  }

} // end namespace ana
