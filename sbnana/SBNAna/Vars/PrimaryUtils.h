/*
Loop over SRTrueParticle vectors and return variables
*/

#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "TVector3.h"
#include <iostream>

using namespace caf;
using namespace std;

namespace ana{

namespace PrimaryUtil{

  // Interaction
  double NeutrinoE(const SRTrueInteractionProxy& true_int);
  int NeutrinoPDG(const SRTrueInteractionProxy& true_int);
  int NeutrinoMode(const SRTrueInteractionProxy& true_int);
  int Target(const SRTrueInteractionProxy& true_int);
  int Npip(const SRTrueInteractionProxy& true_int);
  int Npim(const SRTrueInteractionProxy& true_int);
  int Npi0(const SRTrueInteractionProxy& true_int);

  // Muon
  int MuonIndex(const SRTrueInteractionProxy& true_int);
  double MuonNuCosineTheta(const SRTrueInteractionProxy& true_int);
  double MuonP(const SRTrueInteractionProxy& true_int);
  double MuonPt(const SRTrueInteractionProxy& true_int);
  double MuonCosThBeam(const SRTrueInteractionProxy& true_int);

  // Proton
  int ProtonIndex(const SRTrueInteractionProxy& true_int);
  double ProtonNuCosineTheta(const SRTrueInteractionProxy& true_int);
  double ProtonP(const SRTrueInteractionProxy& true_int);
  double ProtonPt(const SRTrueInteractionProxy& true_int);

  // Muon+Proton
  double CosThMuonProton(const SRTrueInteractionProxy& true_int);

  // TKI
  // https://arxiv.org/abs/1910.08658

  double CalcTKI_deltaPT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTx(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTy(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaalphaT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaphiT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);

  double deltaPT(const SRTrueInteractionProxy& true_int);
  double deltaPTx(const SRTrueInteractionProxy& true_int);
  double deltaPTy(const SRTrueInteractionProxy& true_int);
  double deltaalphaT(const SRTrueInteractionProxy& true_int);
  double deltaphiT(const SRTrueInteractionProxy& true_int);


} // end namespace PrimaryUtil

} // end namespace ana
