#include "sbnana/SBNAna/Cuts/Pi0StudyCuts.h"
#include "sbnana/SBNAna/Vars/Pi0StudyVars.h"

#include <iostream>

namespace ana {

  ////////////////////CRTPMT FlashMatching//////////////////////
  ////////////////////////////////////////////////////////////
  //const SpillCut kIcarus202401CRTPMTVeto([](const caf::SRSpillProxy* spill){
  //  for(const auto& match: spill->crtpmt_matches) {
  //      if(match.flashGateTime > 0 && match.flashGateTime < 1.6 && match.flashClassification == 0)
  //          return true;
  //  }
  //  return false;
  //});
    
  /////////////////Spill Cuts///////////////////////////
  //////////////////////////////////////////////////////

  // Cut on having valid trigger time in approx. beam window pi0study
  const SpillCut kNoSpillCuts([](const caf::SRSpillProxy* sr) {
    return true; // No spill cuts
  });

  const SpillCut kNuMIValidTrigger ( [](const caf::SRSpillProxy *sr) {
    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    return spillTriggerTime > -0.1 && spillTriggerTime < 9.7;
  });

  // Cut on having valid trigger time in approx. beam window pi0study
  const SpillCut kBNBValidTrigger ( [](const caf::SRSpillProxy *sr) {
    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    return spillTriggerTime > -0.1 && spillTriggerTime < 1.6;
  });

  /////////////////Slice Cuts///////////////////////////
  //////////////////////////////////////////////////////
  const Cut kNoSliceCuts ([](const caf::SRSliceProxy* slc) { return true; });

  // reco vertex fiducial volume cut pi0study
  const Cut kNuMIVertexInFV([](const caf::SRSliceProxy* slc) { 
    if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
    return isInFV(slc->vertex.x,slc->vertex.y,slc->vertex.z);
  });

  // reco vertex is contained in active volume cut pi0study
  const Cut kNuMIVertexIsContained([](const caf::SRSliceProxy* slc) {
    if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
    return isContainedVol(slc->vertex.x,slc->vertex.y,slc->vertex.z);
  });

  // reco tagged as clear cosmic cut pi0study
  const Cut kNuMINotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !slc->is_clear_cosmic;
  });

  /////////////Particle Cuts///////////////////////////
  /////////////////////////////////////////////////////
  
  // Muon candidate
  const Cut kNuMIHasMuonCandidate([](const caf::SRSliceProxy* slc) {
    return ( kNuMIMuonCandidateIdx(slc) >= 0 );
  });

  // Proton candidate
  const Cut kNuMIHasProtonCandidate([](const caf::SRSliceProxy* slc) {
    return ( kNuMIProtonCandidateIdx(slc) >= 0 );
  });

  // Dont think I need this cut 
  const Cut kNuMIProtonCandidateRecoPTreshold([](const caf::SRSliceProxy* slc) {
    float p(-5.f);

    if ( kNuMIProtonCandidateIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kNuMIProtonCandidateIdx(slc)).trk;
      p = trk.rangeP.p_proton;
    }

    return (p > 0.0 && p < 1.0);
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
  //

  const Cut kNuMISelection_1muXpi0_Base = kNuMIVertexIsContained
                                          && kNuMIHasMuonCandidate;
                                          //&& kNuMIAllPrimaryHadronsContained;

  const Cut kNuMISelection_1muXpi0_CosmicFilter      = kNuMIVertexIsContained
                                                       //&& kNuMIHasMuonCandidate
                                                       //&& kNuMIAllPrimaryHadronsContained
                                                       && kNuMINotClearCosmic;
  
  const Cut kNuMISelection_1muXpi0_ChargedPionFilter = kNuMIVertexInFV &&
                                                       kNuMINoSecondPrimaryMuonlikeTracks &&
                                                       kNuMIHasMuonCandidate &&
                                                       kNuMIAllPrimaryHadronsContained;

  const Cut kNuMISelection_1muXpi0      = kNuMIVertexInFV &&
                                        kNuMINoSecondPrimaryMuonlikeTracks &&
                                        kNuMINotClearCosmic &&
                                        kNuMIHasMuonCandidate &&
                                        kNuMIAllPrimaryHadronsContained;
  
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

  ////////////////////////////////////////////////////////////////
  // Sidebands
  ////////////////////////////////////////////////////////////////
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

  /// CutType; 1=Signal, 2=pi+- sideband, 3=pi0 sideband, 0=other Gets saved to Tree 
  // From Reco
  const Var kNuMICutType([](const caf::SRSliceProxy* slc) -> double {

    if( kNuMISelection_1muNp0pi(slc) ) return 1; // Full selection
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

  ////////////////////////////////////////////////
  // Signal definitions:
  ////////////////////////////////////////////////
  // Neutrino NC
  const Cut kNuMI_IsSliceNuNC([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false;
    if ( slc->truth.iscc ) return false; // iscc=true -> Slice is not NC
    return true;
  });

  // Not nu matched: i.e. cosmic, or noise, or not well-matched to an interaction
  const Cut kNuMI_IsSlcNotNu([](const caf::SRSliceProxy* slc) {
    return ( slc->truth.index < 0 );
  });

  /// \ref  #1 SIGNAL! Check 1mu1Pi0X using vector of primaries pi0study
  bool Is_1mu_0pipm_1pi0(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
    //From SRTrueInteraction
    if ( true_int.index < 0 ) return false;
    if ( abs(true_int.pdg) != 14
         || !true_int.iscc
         || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
         || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
         || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
         //|| (true_int.time < 0.0 || true_int.time > 9.7) // NuMI beam window
        )
      return false; // not signal

    //From SRTrueParticle from rec.mc.nu.prim
    unsigned int nMu(0), nPi(0), nPi0(0); //, nPhoton(0);
    for ( auto const& prim : true_int.prim ) {
      //if ( prim.start_process != 0 ) continue;
      if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
      //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      double energy = prim.genE;
      int daughters = prim.daughters.size();

      if ( abs(prim.pdg) == 13 && prim.start_process == 0 && momentum > 0.226) nMu+=1;
      if ( abs(prim.pdg) == 211 && prim.start_process == 0 && energy > 0.025) nPi+=1;
      if ( abs(prim.pdg) == 111 && daughters == 2 && prim.start_process == 0) nPi0+=1;
      //if ( abs(prim.pdg) == 22 && prim.start_process == 3 && energy > 0.020) nPhoton+=1;
    }
    if ( nMu==1 && nPi0==1 && nPi==0 ) return true;

    return false;
  }

    /// \ref 1mu_0pi_1pi0 OOPS NOT SIGNAL! Check 1mu1Pi0X using vector of primaries pi0study
    bool Is_1mu_0pipm_1pi0_OOPS(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
      //From SRTrueInteraction
      if ( true_int.index < 0 ) return false;
      if ( abs(true_int.pdg) != 14
            || !true_int.iscc
            || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
            || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
            || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
          )
        return false; // not signal
  
      //From SRTrueParticle from rec.mc.nu.prim
      unsigned int nMu(0), nPi(0), nPi0(0);//, nPhoton(0);
      for ( auto const& prim : true_int.prim ) {
        //if ( prim.start_process != 0 ) continue;
        if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
        //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut
  
        //double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
        //double energy = prim.genE;
        int daughters = prim.daughters.size();
  
        if ( abs(prim.pdg) == 13 && prim.start_process == 0 ) nMu+=1;
        if ( abs(prim.pdg) == 211 && prim.start_process == 0 ) nPi+=1;
        if ( abs(prim.pdg) == 111 && daughters == 2 && prim.start_process == 0) nPi0+=1;
        //if ( abs(prim.pdg) == 22 && prim.start_process == 3 ) nPhoton+=1;
      }
      if ( nMu==1 && nPi0==1 && nPi==0 ) return true;
  
      return false;
    }
    

    /// \ref 1mu_0pi_1pi0 OOFV NOT SIGNAL! Check 1mu1Pi0X using vector of primaries pi0study
  bool Is_1mu_0pipm_1pi0_OOFV(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
    //From SRTrueInteraction
    if ( true_int.index < 0 ) return false;
    if ( abs(true_int.pdg) != 14
        || !true_int.iscc
        || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
        || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
        )
      return false; // not signal

    //From SRTrueParticle from rec.mc.nu.prim
    unsigned int nMu(0), nPi(0), nPi0(0);//, nPhoton(0);
    for ( auto const& prim : true_int.prim ) {
      //if ( prim.start_process != 0 ) continue;
      //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      double energy = prim.genE;
      int daughters = prim.daughters.size();

      if ( abs(prim.pdg) == 13 && prim.start_process == 0 && momentum > 0.226) nMu+=1;
      if ( abs(prim.pdg) == 211 && prim.start_process == 0 && energy > 0.025) nPi+=1;
      if ( abs(prim.pdg) == 111 && daughters == 2 && prim.start_process == 0) nPi0+=1;
      //if ( abs(prim.pdg) == 22 && prim.start_process == 3 && energy > 0.020) nPhoton+=1;
    }
    if ( nMu==1 && nPi0==1 && nPi==0 ) return true;

    return false;
  }

  /// \ref 1mu_Npipm_1pi0 NOT SIGNAL!
  bool Is_1mu_Npipm_1pi0(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
    //From SRTrueInteraction
    if ( true_int.index < 0 ) return false;
    if ( abs(true_int.pdg) != 14
         || !true_int.iscc
         || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
         || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
         || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
        )
      return false; // not signal

    //From SRTrueParticle from rec.mc.nu.prim
    unsigned int nMu(0), nPi(0), nPi0(0);//, nPhoton(0);
    for ( auto const& prim : true_int.prim ) {
      //if ( prim.start_process != 0 ) continue;
      if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
      //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      double energy = prim.genE;
      int daughters = prim.daughters.size();

      if ( abs(prim.pdg) == 13 && prim.start_process == 0 && momentum > 0.226) nMu+=1;
      if ( abs(prim.pdg) == 211 && prim.start_process == 0 && energy > 0.025) nPi+=1;
      if ( abs(prim.pdg) == 111 && daughters == 2 && prim.start_process == 0) nPi0+=1;
      //if ( abs(prim.pdg) == 22 && prim.start_process == 3 && energy > 0.020) nPhoton+=1;
    }
    if ( nMu==1 && nPi >= 1 && nPi0==1 ) return true;

    return false;
  }

  //1mu_Npi_0pi0
  /// \ref 1mu_Npipm_1pi0 NOT SIGNAL!
  bool Is_1mu_Npipm_0pi0(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
    //From SRTrueInteraction
    if ( true_int.index < 0 ) return false;
    if ( abs(true_int.pdg) != 14
        || !true_int.iscc
        || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
        || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
        || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
        )
      return false; // not signal

    //From SRTrueParticle from rec.mc.nu.prim
    unsigned int nMu(0), nPi(0), nPi0(0);//, nPhoton(0);
    for ( auto const& prim : true_int.prim ) {
      //if ( prim.start_process != 0 ) continue;
      if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
      //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      double energy = prim.genE;
      int daughters = prim.daughters.size();

      if ( abs(prim.pdg) == 13 && prim.start_process == 0 && momentum > 0.226) nMu+=1;
      if ( abs(prim.pdg) == 211 && prim.start_process == 0 && energy > 0.025) nPi+=1;
      if ( abs(prim.pdg) == 111 && daughters == 2 && prim.start_process == 0) nPi0+=1;
      //if ( abs(prim.pdg) == 22 && prim.start_process == 3 && energy > 0.020) nPhoton+=1;
    }
    if ( nMu==1 && nPi >= 1 && nPi0==0 ) return true;

    return false;
  }

  /// \ref 1mu_Npi0 NOT SIGNAL!
  bool Is_1mu_Npi0(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
    //From SRTrueInteraction
    if ( true_int.index < 0 ) return false;
    if ( abs(true_int.pdg) != 14
        || !true_int.iscc
        || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
        || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
        || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
        )
      return false; // not signal

    //From SRTrueParticle from rec.mc.nu.prim
    unsigned int nMu(0), nPi(0), nPi0(0);//, nPhoton(0);
    for ( auto const& prim : true_int.prim ) {
      //if ( prim.start_process != 0 ) continue;
      if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
      //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      double energy = prim.genE;
      int daughters = prim.daughters.size();

      if ( abs(prim.pdg) == 13 && prim.start_process == 0 && momentum > 0.226) nMu+=1;
      if ( abs(prim.pdg) == 211 && prim.start_process == 0 && energy > 0.025) nPi+=1;
      if ( abs(prim.pdg) == 111 && daughters == 2 && prim.start_process == 0) nPi0+=1;
      //if ( abs(prim.pdg) == 22 && prim.start_process == 3 && energy > 0.020) nPhoton+=1;
    }
    if ( nMu==1 && nPi0>=1 ) return true;

    return false;
  }

    /// \ref 0mu_1pi0 NOT SIGNAL! IS NC
    bool Is_0mu_1pi0(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
      //From SRTrueInteraction
      if ( true_int.index < 0 ) return false;
      if ( abs(true_int.pdg) != 14
          || true_int.iscc // is NC
          || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
          || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
          || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
          )
        return false; // not signal
  
      //From SRTrueParticle from rec.mc.nu.prim
      unsigned int nMu(0), nPi(0), nPi0(0);//, nPhoton(0);
      for ( auto const& prim : true_int.prim ) {
        //if ( prim.start_process != 0 ) continue;
        if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
        //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut
  
        double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
        double energy = prim.genE;
        int daughters = prim.daughters.size();
  
        if ( abs(prim.pdg) == 13 && prim.start_process == 0 && momentum > 0.226) nMu+=1;
        if ( abs(prim.pdg) == 211 && prim.start_process == 0 && energy > 0.025) nPi+=1;
        if ( abs(prim.pdg) == 111 && daughters == 2 && prim.start_process == 0) nPi0+=1;
        //if ( abs(prim.pdg) == 22 && prim.start_process == 3 && energy > 0.020) nPhoton+=1;
      }
      if ( nMu==0 && nPi0 == 1 ) return true;
  
      return false;
    }


  /// \ref 1mu_Npi_1Pi0 BACKGROUND Check 1muNPi0X using vector of primaries pi0study background
  bool Is1muNPi0X(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
    //From SRTrueInteraction
    if ( true_int.index < 0 ) return false;
    if ( abs(true_int.pdg) != 14
         || !true_int.iscc
         || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
         || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
         || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
        )
      return false; // not signal

    //From SRTrueParticle from rec.mc.nu.prim
    unsigned int nMu(0), nPi(0), nPi0(0);
    for ( auto const& prim : true_int.prim ) {
      //if ( prim.start_process != 0 ) continue;
      
      if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
      //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      double energy = prim.genE;
      int daughters = prim.daughters.size();

      if ( abs(prim.pdg) == 13 && momentum > 0.226 && prim.start_process == 0 ) nMu+=1;
      if ( abs(prim.pdg) == 211 && prim.start_process == 0 && energy > 0.025 ) nPi+=1;
      if ( abs(prim.pdg) == 111 && daughters == 2 ) nPi0+=1;
    }
    if ( nMu==1 && nPi==0 && nPi0 > 1) return true; //1mu_Npi0

    return false;
  }

    /// \ref BACKGROUND Check 1mu0Pi0X using vector of primaries pi0study background
    bool Is1mu0Pi0X(const caf::Proxy<caf::SRTrueInteraction>& true_int){ //Either from rec.slc.truth or rec.mc.nu
      //From SRTrueInteraction
      if ( true_int.index < 0 ) return false;
      if ( abs(true_int.pdg) != 14
           || !true_int.iscc
           || std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z)
           || !isInFV(true_int.position.x, true_int.position.y, true_int.position.z)
           || std::isnan(true_int.prod_vtx.x) || std::isnan(true_int.prod_vtx.y) || std::isnan(true_int.prod_vtx.z)
          )
        return false; // not signal
  
      //From SRTrueParticle from rec.mc.nu.prim
      unsigned int nMu(0), nPi(0), nPi0(0);
      for ( auto const& prim : true_int.prim ) {
        //if ( prim.start_process != 0 ) continue;
        
        if (!isInFV(prim.start.x, prim.start.y, prim.start.z)) continue; // FV cut
        //if (!isInFV(prim.gen.x, prim.gen.y, prim.gen.z)) continue; // FV cut
  
        double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
        //double energy = prim.genE;
        int daughters = prim.daughters.size();
  
        if ( abs(prim.pdg) == 13 && momentum > 0.226 && prim.start_process == 0 ) nMu+=1;
        if ( abs(prim.pdg) == 211 && prim.start_process == 0 ) nPi+=1;
        if ( abs(prim.pdg) == 111 && daughters == 2 ) nPi0+=1;
      }
      if ( nMu==1 && nPi0 == 0 && nPi==0 ) return true;
  
      return false;
    }


  /// \ref Check 1muNp0pi using vector of primaries
  bool Is1muNp0pi(const caf::Proxy<caf::SRTrueInteraction>& true_int, bool ApplyProtonPCut){

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
      if ( abs(prim.pdg) == 13 && momentum > 0.226 ) nMu+=1;

      bool PassProtonPCut = (momentum > 0.4 && momentum < 1.);
      if ( abs(prim.pdg) == 2212 && (ApplyProtonPCut ? PassProtonPCut : true)  ) nP+=1;
      if ( abs(prim.pdg) == 111 || abs(prim.pdg) == 211 ) nPi+=1;
    }
    if ( nMu!=1 || nP==0 || nPi > 0 ) return false;

    return true;

  }

  // CC signal: No kinematic threshold on proton
  const Cut kNuMI_1muNp0piStudy_Signal_NoContainment([](const caf::SRSliceProxy* slc) {
    return Is1muNp0pi(slc->truth, false);
  });

  // CC but not signal: No kinematic threshold on proton
  const Cut kNuMI_1muNp0piStudy_OtherNuCC_NoContainment([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false; // not Nu
    if ( !slc->truth.iscc ) return false; // not CC

    if ( kNuMI_1muNp0piStudy_Signal_NoContainment(slc) ) return false; // covered by signal
    return true;
  });

  // CC signal: With kinematic threshold of 400MeV/c to 1 GeV/c on proton
  const Cut kNuMI_1muNp0piStudy_Signal_NoContainment_ProtonThreshold([](const caf::SRSliceProxy* slc) {
    return Is1muNp0pi(slc->truth, true);
  });

  // OneMu1Pi0 CC SIGNAL #1
  const Cut kNuMI_1mu_0pipm_1pi0_Signal([](const caf::SRSliceProxy* slc) {
    return Is_1mu_0pipm_1pi0(slc->truth); //pi0study
    //return Is1mu1Pi0X(slc->truth); //pi0study OLD 
  });

  // BROKEN OneMu1Pi0 CC SIGNAL #2
  const Cut kNuMI_1mu_0pipm_1pi0_OOPS([](const caf::SRSliceProxy* slc) {
    if ( kNuMI_1mu_0pipm_1pi0_Signal(slc) ) return false; // covered by signal
    return Is_1mu_0pipm_1pi0_OOPS(slc->truth); //pi0study
  });

   // BROKEN OneMu1Pi0 CC SIGNAL #3
  const Cut kNuMI_1mu_0pipm_1pi0_OOFV([](const caf::SRSliceProxy* slc) {
    if ( kNuMI_1mu_0pipm_1pi0_Signal(slc) ) return false; // covered by signal
    return Is_1mu_0pipm_1pi0_OOFV(slc->truth); //pi0study
  });

  // Background Is_1mu_Npipm_1pi0 #4
  const Cut kNuMI_Is_1mu_Npipm_1pi0([](const caf::SRSliceProxy* slc) {
    if ( kNuMI_1mu_0pipm_1pi0_Signal(slc) ) return false; // covered by signal
    return Is_1mu_Npipm_1pi0(slc->truth);
  });

   // Background Is_1mu_Npipm_0pi0 #5
   const Cut kNuMI_Is_1mu_Npipm_0pi0([](const caf::SRSliceProxy* slc) {
    if ( kNuMI_1mu_0pipm_1pi0_Signal(slc) ) return false; // covered by signal
    return Is_1mu_Npipm_0pi0(slc->truth);
  });

    // Background Is_1mu_Npi0 #6
    const Cut kNuMI_Is_1mu_Npi0([](const caf::SRSliceProxy* slc) {
      if ( kNuMI_1mu_0pipm_1pi0_Signal(slc) ) return false; // covered by signal
      return Is_1mu_Npi0(slc->truth);
    });

    // Background Is_0mu_Npi0 #7
    const Cut kNuMI_Is_0mu_1pi0([](const caf::SRSliceProxy* slc) {
      if ( kNuMI_1mu_0pipm_1pi0_Signal(slc) ) return false; // covered by signal
      return Is_0mu_1pi0(slc->truth);
    });

    //Background Is_Other_nu #8
    const Cut kNuMI_Is_Other_nu([](const caf::SRSliceProxy* slc) {
      if ( slc->truth.index < 0 ) return false; // not Nu
      return true;
    });

  // CC OneMu0Pi0. THIS IS NOT SIGNAL!!
  //const Cut kNuMI_1mu1Pi0XStudy_CC1mu0pi0X([](const caf::SRSliceProxy* slc) {
  //  if ( kNuMI_1mu1Pi0XStudy_Signal(slc) ) return false; // covered by signal
  //  return Is1mu0Pi0X(slc->truth);
  //});

  // CC but not signal: With kinematic threshold of 400MeV/c to 1 GeV/c on proton
  /*
  const Cut kNuMI_1mu1Pi0XStudy_OtherCCNuMu([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false; // not Nu
    if ( !slc->truth.iscc ) return false; // not CC
    if ( !(abs(slc->truth.pdg) == 14) ) return false; // not NuMu
    if ( kNuMI_1mu1Pi0XStudy_Signal(slc) ) return false; // covered by signal
    if ( kNuMI_1mu1Pi0XStudy_CC1muNpi0X(slc) ) return false; // covered by prev background
    if ( kNuMI_1mu1Pi0XStudy_CC1mu0pi0X(slc) ) return false; // covered by prev background
    return true;
  });
    */

  // CC but its a nue
  const Cut kNuMI_1mu1Pi0XStudy_CCNue([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false; // not Nu
    if ( !slc->truth.iscc ) return false; // not CC
    if ( !(abs(slc->truth.pdg) == 12) ) return false; // not NuE
    return true;
  });

  // NC
  const Cut kNuMI_1mu1Pi0XStudy_NC([](const caf::SRSliceProxy* slc) {
    return kNuMI_IsSliceNuNC(slc);
  });

  // Not nu matched: i.e. cosmic, or noise, or not well-matched to an interaction
  const Cut kNuMI_1mu1Pi0XStudy_NotNu([](const caf::SRSliceProxy* slc) {
    return kNuMI_IsSlcNotNu(slc);
  });

  /// CutType; 
  // 1=Signal,
  // 2=1mu1pi0x_OOPS,
  // 3=1mu1pi0x_OOFV,
  // 4=CC1muNpi0X with N > 1 No PiM or PiP, 
  // 5=CC1mu0pi0X, 
  // 6=OtherNuMuCC, 
  // 7=CCNue, 
  // 8=NC, 
  // 9=NotNu, 
  /// can be expanded further
  const Var kNuMISliceSignalType([](const caf::SRSliceProxy* slc) -> int {
    if ( kNuMI_1mu_0pipm_1pi0_Signal(slc) ) return 1; //Signal definition from truth
    else if ( kNuMI_1mu_0pipm_1pi0_OOPS(slc) ) return 2; //Signal definition from truth
    else if ( kNuMI_1mu_0pipm_1pi0_OOFV(slc) ) return 3; //Signal definition from truth
    else if ( kNuMI_Is_1mu_Npipm_1pi0(slc) ) return 4;
    else if ( kNuMI_Is_1mu_Npipm_0pi0(slc) ) return 5;
    else if ( kNuMI_Is_1mu_Npi0(slc) ) return 6;
    else if ( kNuMI_Is_0mu_1pi0(slc) ) return 7;
    else if ( kNuMI_Is_Other_nu(slc) ) return 8;
    else return 0;
  });

}
