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
  extern const Var MuonTrackP;
  extern const Var MuonTrackNuMICosineTheta;
  extern const Var MuonTrackNuMIToVtxCosineTheta;
  // truth match
  extern const Var MuonTrackTruthLength;
  extern const Var MuonTrackTruthP;
  extern const Var MuonTrackTruthNuMICosineTheta;

  extern const Var ProtonTrackIndex;
  extern const Var ProtonTrackLength;
  extern const Var ProtonTrackP;
  extern const Var ProtonTrackNuMICosineTheta;
  extern const Var ProtonTrackNuMIToVtxCosineTheta;
  // truth match
  extern const Var ProtonTrackTruthLength;
  extern const Var ProtonTrackTruthP;
  extern const Var ProtonTrackTruthNuMICosineTheta;

  extern const Var MuonProtonCosineTheta;

  namespace Aux{
    extern const Var RelaxedMuonTrackIndex;
    extern const Var RelaxedMuonTrackLength;
    extern const Var RelaxedMuonTrackP;
    extern const Var RelaxedMuonTrackDirX;
    extern const Var RelaxedMuonTrackDirY;
    extern const Var RelaxedMuonTrackDirZ;
    extern const Var RelaxedMuonTrackTrackScore;
    extern const Var RelaxedMuonTrackChi2Muon;
    extern const Var RelaxedMuonTrackChi2Proton;
    extern const Var RelaxedMuonTrackChi2MuonCollection;
    extern const Var RelaxedMuonTrackChi2ProtonCollection;
    extern const Var RelaxedMuonTrackCustomChi2MuonCollection;
    extern const Var RelaxedMuonTrackCustomChi2ProtonCollection;
    extern const MultiVar RelaxedMuonTrackCollectionRR;
    extern const MultiVar RelaxedMuonTrackCollectiondEdX;
    extern const MultiVar RelaxedMuonTrackCollectiondQdX;
    // truth match
    extern const Var RelaxedMuonTrackTruthLength;
    extern const Var RelaxedMuonTrackTruthStartProcess;
    extern const Var RelaxedMuonTrackTruthEndProcess;

    extern const Var RelaxedProtonTrackIndex;
    extern const Var RelaxedProtonTrackLength;
    extern const Var RelaxedProtonTrackDirX;
    extern const Var RelaxedProtonTrackDirY;
    extern const Var RelaxedProtonTrackDirZ;
    extern const Var RelaxedProtonTrackTrackScore;
    extern const Var RelaxedProtonTrackChi2Muon;
    extern const Var RelaxedProtonTrackChi2Proton;
    extern const Var RelaxedProtonTrackChi2MuonCollection;
    extern const Var RelaxedProtonTrackChi2ProtonCollection;
    extern const Var RelaxedProtonTrackCustomChi2MuonCollection;
    extern const Var RelaxedProtonTrackCustomChi2ProtonCollection;
    extern const Var RelaxedProtonTrackNHitCollection;
    extern const MultiVar RelaxedProtonTrackCollectionRR;
    extern const MultiVar RelaxedProtonTrackCollectiondEdX;
    extern const MultiVar RelaxedProtonTrackCollectiondQdX;
    // truth match
    extern const Var RelaxedProtonTrackTruthLength;
    extern const Var RelaxedProtonTrackTruthStartProcess;
    extern const Var RelaxedProtonTrackTruthEndProcess;

  }

} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
