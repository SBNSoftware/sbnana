/*
Loop over SRTrueParticle vectors and return variables
*/

#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "TVector3.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ana{

namespace PrimaryUtil{

  using Primary = caf::Proxy<std::vector<caf::SRTrueParticle>>;
  using TrueInteraction = caf::Proxy<caf::SRTrueInteraction>;

  // Interaction
  double NeutrinoE(const TrueInteraction& true_int);
  double NeutrinoPDG(const TrueInteraction& true_int);
  double NeutrinoMode(const TrueInteraction& true_int);
  double Target(const TrueInteraction& true_int);

  // Muon
  int MuonIndex(const TrueInteraction& true_int);
  double MuonNuCosineTheta(const TrueInteraction& true_int);
  double MuonP(const TrueInteraction& true_int);
  double MuonPt(const TrueInteraction& true_int);
  double MuonCosThBeam(const TrueInteraction& true_int);

  // Proton
  int ProtonIndex(const TrueInteraction& true_int);
  double ProtonNuCosineTheta(const TrueInteraction& true_int);
  double ProtonP(const TrueInteraction& true_int);
  double ProtonPt(const TrueInteraction& true_int);

  // Muon+Proton
  double CosThMuonProton(const TrueInteraction& true_int);

  // TKI
  // https://arxiv.org/abs/1910.08658

  double CalcTKI_deltaPT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTx(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTy(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);

  double deltaPT(const TrueInteraction& true_int);
  double deltaPTx(const TrueInteraction& true_int);
  double deltaPTy(const TrueInteraction& true_int);


} // end namespace PrimaryUtil

} // end namespace ana
