/*
Loop over SRTrueParticle vectors and return variables
*/

#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "TVector3.h"
#include <iostream>

#define M_MUON 0.1057
#define M_CHARGEDPION 0.13957039
#define M_NEUTRALPION 0.1349768
#define M_PIZERO 0.1349768
#define M_PROTON 0.938272
#define M_NEUTRON 0.939565
#define M_ELECTRON 0.00051
#define E_EffNuclB 0.040

namespace ana{

namespace PrimaryUtil{

  // Interaction
  double NeutrinoE_True(const caf::SRTrueInteractionProxy& true_int);
  int NeutrinoPDG_True(const caf::SRTrueInteractionProxy& true_int);
  int NeutrinoMode_True(const caf::SRTrueInteractionProxy& true_int);
  int IsCC_True(const caf::SRTrueInteractionProxy& true_int);
  int Target_True(const caf::SRTrueInteractionProxy& true_int);
  int NProton_True(const caf::SRTrueInteractionProxy& true_int);
  int NNeutron_True(const caf::SRTrueInteractionProxy& true_int);
  int Npip_True(const caf::SRTrueInteractionProxy& true_int);
  int Npip_True_Any(const caf::SRTrueInteractionProxy& true_int);
  int Npim_True(const caf::SRTrueInteractionProxy& true_int);
  int Npim_True_Any(const caf::SRTrueInteractionProxy& true_int);
  int Npi0_True(const caf::SRTrueInteractionProxy& true_int);
  int Npi0_True_Any(const caf::SRTrueInteractionProxy& true_int);
  double Q2_True(const caf::SRTrueInteractionProxy& true_int);
  double q0_True(const caf::SRTrueInteractionProxy& true_int);
  double q3_True(const caf::SRTrueInteractionProxy& true_int);
  double w_True(const caf::SRTrueInteractionProxy& true_int);

  // Muon
  int MuonIndex_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonNuCosineTheta_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonCosThBeam_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonP_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonPt_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonKE_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonLength_True(const caf::SRTrueInteractionProxy& true_int);
  int MuonContained_True(const caf::SRTrueInteractionProxy& true_int); // 1: contained, 0: not contained (-1: muon not found)

  // (Leading) Proton
  int ProtonIndex_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonNuCosineTheta_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonCosThBeam_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonP_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonPt_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonKE_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonLength_True(const caf::SRTrueInteractionProxy& true_int);

  // Muon+Proton
  double CosThMuonProton_True(const caf::SRTrueInteractionProxy& true_int);

  // TKI
  // https://arxiv.org/abs/1910.08658

  double CalcTKI_deltaPT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTx(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTy(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaalphaT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaphiT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);

  double deltaPT_True(const caf::SRTrueInteractionProxy& true_int);
  double deltaPTx_True(const caf::SRTrueInteractionProxy& true_int);
  double deltaPTy_True(const caf::SRTrueInteractionProxy& true_int);
  double deltaalphaT_True(const caf::SRTrueInteractionProxy& true_int);
  double deltaphiT_True(const caf::SRTrueInteractionProxy& true_int);

  // (Leading) Charged pion
  int ChargedPionIndex_True(const caf::SRTrueInteractionProxy& true_int);
  double ChargedPionKE_True(const caf::SRTrueInteractionProxy& true_int);


} // end namespace PrimaryUtil

} // end namespace ana
