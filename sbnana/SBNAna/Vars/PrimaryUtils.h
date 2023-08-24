/*
Loop over SRTrueParticle vectors and return variables
*/

#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "TVector3.h"
#include <iostream>

namespace ana{

namespace PrimaryUtil{

  // Interaction
  double NeutrinoE_True(const caf::SRTrueInteractionProxy& true_int);
  int NeutrinoPDG_True(const caf::SRTrueInteractionProxy& true_int);
  int NeutrinoMode_True(const caf::SRTrueInteractionProxy& true_int);
  int Target_True(const caf::SRTrueInteractionProxy& true_int);
  int Npip_True(const caf::SRTrueInteractionProxy& true_int);
  int Npim_True(const caf::SRTrueInteractionProxy& true_int);
  int Npi0_True(const caf::SRTrueInteractionProxy& true_int);

  // Muon
  int MuonIndex_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonNuCosineTheta_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonP_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonPt_True(const caf::SRTrueInteractionProxy& true_int);
  double MuonCosThBeam_True(const caf::SRTrueInteractionProxy& true_int);

  // Proton
  int ProtonIndex_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonNuCosineTheta_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonP_True(const caf::SRTrueInteractionProxy& true_int);
  double ProtonPt_True(const caf::SRTrueInteractionProxy& true_int);

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


} // end namespace PrimaryUtil

} // end namespace ana
