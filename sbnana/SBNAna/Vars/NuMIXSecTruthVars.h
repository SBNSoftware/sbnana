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

  // Neutrino/interaction

  const TruthVar kTruth_NeutrinoE = SIMPLETRUTHVAR(E);
  const TruthVar kTruth_NeutrinoPDG = SIMPLETRUTHVAR(pdg);
  const TruthVar kTruth_NeutrinoMode = SIMPLETRUTHVAR(genie_mode);
  const TruthVar kTruth_IsCC = SIMPLETRUTHVAR(iscc);
  const TruthVar kTruth_Target = SIMPLETRUTHVAR(targetPDG);
  const TruthVar kTruth_Q2 = SIMPLETRUTHVAR(Q2);
  const TruthVar kTruth_q0 = SIMPLETRUTHVAR(q0_lab);
  const TruthVar kTruth_q3 = SIMPLETRUTHVAR(modq_lab);
  const TruthVar kTruth_w = SIMPLETRUTHVAR(w);
  extern const TruthVar kTruth_NProton_Primary;
  extern const TruthVar kTruth_NNeutron_Primary;
  extern const TruthVar kTruth_Npip_Primary;
  extern const TruthVar kTruth_Npip_All;
  extern const TruthVar kTruth_Npim_Primary;
  extern const TruthVar kTruth_Npim_All;
  extern const TruthVar kTruth_Npi0_Primary;
  extern const TruthVar kTruth_Npi0_All;
  extern const TruthVar kTruth_IsFHC; // 0: RHC, 1: FHC; note) dummy for now. always return 1

  // Muon

  extern const TruthVar kTruth_MuonIndex;
  extern const TruthVar kTruth_MuonNuCosineTheta;
  extern const TruthVar kTruth_MuonCosThBeam;
  extern const TruthVar kTruth_MuonP;
  extern const TruthVar kTruth_MuonPt;
  extern const TruthVar kTruth_MuonKE;
  extern const TruthVar kTruth_MuonLength;
  extern const TruthVar kTruth_MuonContained; // 1: contained, 0: not contained (-1: muon not found)

  // (Leading) Proton

  extern const TruthVar kTruth_ProtonIndex;
  extern const TruthVar kTruth_ProtonNuCosineTheta;
  extern const TruthVar kTruth_ProtonCosThBeam;
  extern const TruthVar kTruth_ProtonP;
  extern const TruthVar kTruth_ProtonPt;
  extern const TruthVar kTruth_ProtonKE;
  extern const TruthVar kTruth_ProtonLength;

  // Muon+Proton
  extern const TruthVar kTruth_CosThMuonProton;

  // TKI calculator using momentum vectors
  // https://arxiv.org/abs/1910.08658
  // Can be used for both reco and truth
  double CalcTKI_deltaPT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTx(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaPTy(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaalphaT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  double CalcTKI_deltaphiT(const TVector3 vec_p_mu, const TVector3 vec_p_pro, const TVector3 vec_p_nu);
  // TKI
  extern const TruthVar kTruth_deltaPT;
  extern const TruthVar kTruth_deltaPTx;
  extern const TruthVar kTruth_deltaPTy;
  extern const TruthVar kTruth_deltaalphaT;
  extern const TruthVar kTruth_deltaphiT;

  // (Leading) Charged pion

  extern const TruthVar kTruth_ChargedPionIndex;
  extern const TruthVar kTruth_ChargedPionKE;

} // end namespace ana
