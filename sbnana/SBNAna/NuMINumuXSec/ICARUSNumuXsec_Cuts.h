#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  // SpillCut

  // - CRTPMT matching

  extern const SpillCut spillcutHasValidFlash;
  extern const SpillCut spillcutHasInTimeFlash;
  extern const SpillCut spillcutHasOneTimeFlash;
  extern const SpillCut spillcutAllEnteringByCRTPMT;
  extern const SpillCut spillcutHasCRTPMTDirtByID;

  // - Spill nu
  extern const SpillCut spillcutHasTrueMuContained;
  extern const SpillCut spillcutHasTrueMuExiting;
  // - Spill with CC neutrino
  extern const SpillCut spillcutHasCCNeutrino;
  extern const SpillCut spillcutVisECutGT200MeV;
  extern const SpillCut spillcutVisECutGT400MeV;
  extern const SpillCut spillcutVisECutGT600MeV;

  // Cut (slice)

  // - Interaction
  //   - Neutrino flavor
  extern const Cut cutIsNuMu;
  extern const Cut cutIsNuE;
  //   - CC vs NC
  extern const Cut cutIsCC;
  extern const Cut cutIsNC;
  //   - GENIE Interaction code
  extern const Cut cutIsQE;
  extern const Cut cutIsRes;
  extern const Cut cutIsDIS;
  extern const Cut cutIsCoh;
  extern const Cut cutIsCohElastic;
  extern const Cut cutIsElectronScattering;
  extern const Cut cutIsIMDAnnihilation;
  extern const Cut cutIsInverseBetaDecay;
  extern const Cut cutIsGlashowResonance;
  extern const Cut cutIsAMNuGamma;
  extern const Cut cutIsMEC;
  extern const Cut cutIsDiffractive;
  extern const Cut cutIsEM;
  extern const Cut cutIsWeakMix;
  extern const Cut cutIsUnknownInteractionType1;
  extern const Cut cutIsUnknownInteractionType2;
  extern const Cut cutIsUnknownInteractionType3;
  //   - NuMu-CC categories
  extern const Cut cutIsNuMuCC;
  extern const Cut cutIsNuMuCCQE;
  extern const Cut cutIsNuMuCCRes;
  extern const Cut cutIsNuMuCCMEC;
  extern const Cut cutIsNuMuCCDIS;
  extern const Cut cutIsNuMuCCCoh;
  extern const Cut cutIsNuMuCCCohElastic;
  //   - NuMu-NC categories
  extern const Cut cutIsNuMuNC;
  //   - NuE
  extern const Cut cutIsNuECC;

  // - FV

  extern const Cut cutRecoVertexFV;
  extern const Cut cutTruthVertexActiveVolume;
  extern const Cut cutRFiducial;

  // - Pandora decision

  extern const Cut cutNotClearCosmic;
  extern const Cut cutCRLongestTrackDirY;
  extern const Cut cutCRLongestTrackDirYHard;

  // - Flash matching

  extern const Cut cutFMScore;
  extern const Cut cutFMTime;

}
