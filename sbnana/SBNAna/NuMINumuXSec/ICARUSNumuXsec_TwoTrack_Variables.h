#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  extern const Var MuonTrackIndex;
  extern const Var MuonTrackLength;
  extern const Var MuonTrackLengthMatchMuon;
  extern const Var MuonTrackLengthMatchPionPlus;
  extern const Var MuonTrackLengthMatchPionMinus;
  extern const Var MuonTrackLengthMatchProton;
  extern const Var MuonTrackLengthMatchElse;
  extern const Var MuonTrackP;
  extern const Var MuonTrackPMatchMuon;
  extern const Var MuonTrackPMatchPionPlus;
  extern const Var MuonTrackPMatchPionMinus;
  extern const Var MuonTrackPMatchProton;
  extern const Var MuonTrackPMatchElse;
  extern const Var MuonTrackNuMICosineTheta;
  extern const Var MuonTrackNuMIToVtxCosineTheta;

  extern const Var ProtonTrackIndex;
  extern const Var ProtonTrackP;
  extern const Var ProtonTrackPMatchMuon;
  extern const Var ProtonTrackPMatchPionPlus;
  extern const Var ProtonTrackPMatchPionMinus;
  extern const Var ProtonTrackPMatchProton;
  extern const Var ProtonTrackPMatchElse;
  extern const Var ProtonTrackNuMICosineTheta;
  extern const Var ProtonTrackNuMIToVtxCosineTheta;

  namespace Aux{
    extern const Var RelaxedMuonTrackIndex;
    extern const Var RelaxedMuonTrackTrackScore;
    extern const Var RelaxedMuonTrackChi2Muon;
    extern const Var RelaxedMuonTrackChi2Proton;

    extern const Var RelaxedProtonTrackIndex;
    extern const Var RelaxedProtonTrackTrackScore;
    extern const Var RelaxedProtonTrackChi2Muon;
    extern const Var RelaxedProtonTrackChi2Proton;
  }

} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
