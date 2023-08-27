#pragma once

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/PrimaryUtils.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

namespace ana{

  // True
  // - Interaction
  extern const Var kNuMITruePDG; //!< Neutrino pdg
  extern const Var kNuMITrueTarget; //!< Target pdg
  extern const Var kNuMITrueMode; //!< GENIE interaction code (https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360)
  extern const Var kNuMITrueIsCC; //!< IsCC (0:NC, 1:CC, -1:Not neutrino)
  extern const Var kNuMITrueNProton; //!< Number of primary proton
  extern const Var kNuMITrueNNeutron; //!< Number of primary neutron
  extern const Var kNuMITrueNpip; //!< Number of primary pi+
  extern const Var kNuMITrueNpim; //!< Number of primary pi-
  extern const Var kNuMITrueNpi0; //!< Number of primary pi0
  extern const Var kNuMITrueNuE; //!< Neutrino energy
  extern const Var kNuMITrueQ2; //!< Q2
  extern const Var kNuMITrueq0; //!< q0; energy transfer
  extern const Var kNuMITrueq3; //!< q3; momentum transfer
  extern const Var kNuMITruew; //!< w; hadronic mass
  // - Muon
  extern const Var kNuMIMuonTrueKE; //!< True muon kinetic energy
  extern const Var kNuMIMuonNuCosineTheta; //!< True muon cosine angle w.r.t. neutrino
  // - Proton
  extern const Var kNuMIProtonTrueKE; //!< True proton kinetic energy
  extern const Var kNuMIProtonNuCosineTheta; //!< True muon cosine angle w.r.t. neutrino
  // - Signal def
  extern const Var kNuMIIsSignal; //!< 0: not signal, 1: signal

  // Reco
  // - Muon
  extern const Var kNuMIRecoMuonContained; //!< 0: Muon candidate track exiting, 1: Muon candidate track contained (-1: no muon candidate)


}
