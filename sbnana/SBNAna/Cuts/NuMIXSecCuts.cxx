#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"

#include <iostream>

namespace ana {

  // Cut on having valid trigger time in approx. beam window
  const SpillCut kNuMIValidTrigger ( [](const caf::SRSpillProxy *sr) {
    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    return spillTriggerTime > -0.1 && spillTriggerTime < 9.7;
  });

  // reco vertex fiducial volume cut
  const Cut kNuMIVertexInFV([](const caf::SRSliceProxy* slc) {
    if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
    return isInFV(slc->vertex.x,slc->vertex.y,slc->vertex.z);
  });

  // reco tagged as clear cosmic cut
  const Cut kNuMINotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !slc->is_clear_cosmic;
  });


  // Muon candidate
  const Cut kNuMIHasMuonCandidate([](const caf::SRSliceProxy* slc) {
    return ( kNuMIMuonCandidateIdx(slc) >= 0 );
  });

  // Proton candidate
  const Cut kNuMIHasProtonCandidate([](const caf::SRSliceProxy* slc) {
    return ( kNuMIProtonCandidateIdx(slc) >= 0 );
  });

  const Cut kNuMIProtonCandidateRecoPTreshold([](const caf::SRSliceProxy* slc) {
    float p(-5.f);

    if ( kNuMIProtonCandidateIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kNuMIProtonCandidateIdx(slc)).trk;
      p = trk.rangeP.p_proton;
    }

    return (p > 0.4 && p < 1.0);
  });

  // Set of two cuts largely meant to cut charged pions
  const Cut kNuMIAllPrimaryHadronsContained([](const caf::SRSliceProxy* slc) {
    // Considers all track fits, not just track-like PFPs
    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

      if ( !isPrimCandidate ) continue;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if ( !Contained ) return false;
    }

    return true;
  });

  const Cut kNuMINoSecondPrimaryMuonlikeTracks([](const caf::SRSliceProxy* slc) {
    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    std::vector<double> chargedpion_indices = kNuMIChargedPionCandidateIdxs(slc);
    if(chargedpion_indices.size()==0) return true;
    else return false;

  });

  // Cut on showers aiming at rejecting pi0
  const Cut kNuMICutPhotons([](const caf::SRSliceProxy* slc) {

    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    int primaryProtonInd = kNuMIProtonCandidateIdx(slc);
    if ( primaryProtonInd < 0 ) return false;

    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()==0) return true;
    else return false;

  });

  // Base selection common to side-bands (cuts back on number of entries one needs to carry):
  // For now the sidebands revolve around chg pi and pi0
  const Cut kNuMISelection_1muNp_Base = kNuMIVertexInFV && kNuMINotClearCosmic &&
                                        kNuMIHasMuonCandidate && kNuMIHasProtonCandidate &&
                                        kNuMIProtonCandidateRecoPTreshold &&
                                        kNuMIAllPrimaryHadronsContained;

  // Full selection: 1muNp0pi without caring about containment
  const Cut kNuMISelection_1muNp0pi_WithoutShowerCut = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                                       kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                                       kNuMIAllPrimaryHadronsContained && kNuMINoSecondPrimaryMuonlikeTracks;                   /*Esp. chg pi rejection*/

  // Full selection: 1muNp0pi without caring about containment
  const Cut kNuMISelection_1muNp0pi = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                      kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                      kNuMIAllPrimaryHadronsContained && kNuMINoSecondPrimaryMuonlikeTracks &&                 /*Esp. chg pi rejection*/
                                      kNuMICutPhotons;                                                                         /*Esp. pi0 rejection*/

  // Muon containment
  const Cut kNuMIMuonCandidateContained([](const caf::SRSliceProxy* slc) {
    int muonCandidate = kNuMIMuonCandidateIdx(slc);

    if ( muonCandidate < 0 ) return false;
    else {
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      return Contained;
    }

    return false;
  });

  // Reco quality (meant to help ID split tracks)
  const Cut kNuMIRejectSplitMuons([](const caf::SRSliceProxy* slc) {
    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isSplitCandidate = (Atslc > 10.);

      if ( !isSplitCandidate || trk.calo[2].nhit < 5 ) continue;
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( !(Chi2Proton > 60. && Chi2Muon < 30.) ) continue;

      if ( trk.len > 10. ) return false;
    }

    return true;
  });

  const Cut kNuMIHasTwoPhotons([](const caf::SRSliceProxy* slc) {

    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    int primaryProtonInd = kNuMIProtonCandidateIdx(slc);
    if ( primaryProtonInd < 0 ) return false;

    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()>=2) return true;
    else return false;

  });

  /// Pion sideband
  const Cut kNuMIChargedPionSideBand = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                       kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                       kNuMIAllPrimaryHadronsContained && /*Hadron contained*/
                                       !kNuMINoSecondPrimaryMuonlikeTracks && /*Charged pion*/
                                       kNuMICutPhotons; /*Esp. pi0 rejection*/

  const Cut kNuMINeutralPionSideBand = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                       kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                       kNuMIAllPrimaryHadronsContained && /*Hadron contained*/
                                       kNuMINoSecondPrimaryMuonlikeTracks && /*Esp. chg pi rejection*/
                                       !kNuMICutPhotons; /*Neutral pion*/

  /// Sideband with the 2ph cut --> this is a subset of the NeutralPionSideband
  const Cut kNuMINeutralPion2phSideBand = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                          kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                          kNuMIAllPrimaryHadronsContained && /*Hadron contained*/
                                          kNuMINoSecondPrimaryMuonlikeTracks && /*Esp. chg pi rejection*/
                                          kNuMIHasTwoPhotons; /*Neutral pion*/

  /// CutType; 1=Signal, 2=pi+- sideband, 3=pi0 sideband, 0=other
  const Var kNuMICutType([](const caf::SRSliceProxy* slc) -> double {

    if( kNuMISelection_1muNp0pi(slc) ) return 1;
    else if( kNuMIChargedPionSideBand(slc) ) return 2;
    else if( kNuMINeutralPion2phSideBand(slc) ) return 3; // this should be a subset of 4! so if you want 4, use "3 || 4"
    else if( kNuMINeutralPionSideBand(slc) ) return 4;
    else return 0;

  });

  /// CutType without the showers cut in signal selection; 1=Signal, 2=pi+- sideband, 3=pi0 sideband, 0=other
  const Var kNuMICutTypeWithoutShowerCut([](const caf::SRSliceProxy* slc) -> double {

    if( kNuMISelection_1muNp0pi_WithoutShowerCut(slc) ) return 1;
    else if( kNuMIChargedPionSideBand(slc) ) return 2;
    else if( kNuMINeutralPion2phSideBand(slc) ) return 3; // this should be a subset of 4! so if you want 4, use "3 || 4"
    else if( kNuMINeutralPionSideBand(slc) ) return 4;
    else return 0;

  });

  // Signal definitions:
  // Neutrino NC
  const Cut kNuMI_IsSliceNuNC([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false;

    if ( slc->truth.iscc )
      return false; // IS CC

    return true;
  });

  // Not nu matched: i.e. cosmic, or noise, or not well-matched to an interaction
  const Cut kNuMI_IsSlcNotNu([](const caf::SRSliceProxy* slc) {
    return ( slc->truth.index < 0 );
  });
  /// \ref Check 1muNp0pi using vector of primaries
  bool Is1muNp0pi(const caf::Proxy<caf::SRTrueInteraction>& true_int, bool ApplyPhaseSpcaeCut){

    if ( true_int.index < 0 ) return false;

    if ( abs(true_int.pdg) != 14 ||
         !true_int.iscc ||
         std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z) ||
         !isInFV(true_int.position.x, true_int.position.y, true_int.position.z) )
      return false; // not signal

    unsigned int nMu(0), nP(0), nPi(0);
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );

      bool PassMuonPCut = (momentum > 0.226);
      if ( abs(prim.pdg) == 13 && (ApplyPhaseSpcaeCut ? PassMuonPCut : true) ) nMu+=1;

      bool PassProtonPCut = (momentum > 0.4 && momentum < 1.);
      if ( abs(prim.pdg) == 2212 && (ApplyPhaseSpcaeCut ? PassProtonPCut : true)  ) nP+=1;

      if ( abs(prim.pdg) == 111 || abs(prim.pdg) == 211 ) nPi+=1;
    }
    if ( nMu!=1 || nP==0 || nPi > 0 ) return false;

    return true;

  }
  /// \ref Signal without phase space cut
  const Cut kNuMI_1muNp0piStudy_Signal_WithoutPhaseSpaceCut([](const caf::SRSliceProxy* slc) {
    return Is1muNp0pi(slc->truth, false);
  });
  /// \ref Signal with phase space cut = "Signal"
  const Cut kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut([](const caf::SRSliceProxy* slc) {
    return Is1muNp0pi(slc->truth, true);
  });
  const TruthCut kTruthCut_IsSignal([](const caf::SRTrueInteractionProxy* nu) {
    return Is1muNp0pi(*nu, true);
  });
  /// \ref Signal but fails phase space cut = "out of phase space" (OOPS)
  const Cut kNuMI_1muNp0piStudy_Signal_FailPhaseSpaceCut = kNuMI_1muNp0piStudy_Signal_WithoutPhaseSpaceCut && // signal,
                                                           !kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut; // but fail phase cut
  /// \ref CC but NOT "Signal" or "OOPS"
  const Cut kNuMI_1muNp0piStudy_OtherNuCC([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false; // not Nu
    if ( !slc->truth.iscc ) return false; // not CC

    if ( kNuMI_1muNp0piStudy_Signal_WithoutPhaseSpaceCut(slc) ) return false; // covered by signal without phase space cut
    return true;
  });

  /// CutType; 1=Signal, 2=OOPS, 3=OtherCC, 4=NuNC, 5=NotNu
  /// can be expanded further
  const Var kNuMISliceSignalType([](const caf::SRSliceProxy* slc) -> int {
    if ( kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return 1; // Signal (with phase space cut)
    else if ( kNuMI_1muNp0piStudy_Signal_FailPhaseSpaceCut(slc) ) return 2; // Signal but out of phase space cut (OOPS)
    else if ( kNuMI_1muNp0piStudy_OtherNuCC(slc) ) return 3; // CC but not (signal without phase space cut)
    else if ( kNuMI_IsSliceNuNC(slc) ) return 4; // NC
    else if ( kNuMI_IsSlcNotNu(slc) ) return 5; // Not nu-slice Cosmic
    else return 0;
  });

}
