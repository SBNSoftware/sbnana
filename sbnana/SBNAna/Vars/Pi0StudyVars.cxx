#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnana/SBNAna/Vars/Pi0StudyVars.h"
#include "sbnana/SBNAna/Cuts/Pi0StudyCuts.h"

#include "TVector3.h"

namespace ana {

  // Utility functions
  bool isInFV (double x, double y, double z)
  {
    if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
    if ( std::isinf(x) || std::isinf(y) || std::isinf(z) ) return false;

    return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
              ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
            ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
              ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
  }

  bool isContainedVol (double x, double y, double z)
  {
    if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;
    if ( std::isinf(x) || std::isinf(y) || std::isinf(z) ) return false;

    return (( ( x < -61.94 - 5. && x > -358.49 + 5. ) ||
              ( x >  61.94 + 5. && x <  358.49 - 5. )) &&
            ( ( y > -181.86 + 5. && y < 134.96 - 5. ) &&
              ( z > -894.95 + 5. && z < 894.95 - 5. ) ));
  }

  bool IsValidTrkIdx( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
    return slice->reco.npfp > idxTrk;
  }

  bool IsTracklikeTrack( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
    return (!std::isnan(slice->reco.pfp.at(idxTrk).trackScore) && slice->reco.pfp.at(idxTrk).trackScore > 0.45);
  }

  bool IsShowerlike( const caf::SRSliceProxy* slice, const unsigned int idxShw ) {
    return (
            !std::isnan(slice->reco.pfp.at(idxShw).trackScore) && 
            slice->reco.pfp.at(idxShw).trackScore > 0. && 
            slice->reco.pfp.at(idxShw).trackScore <= 0.45
            );
  }

  bool IsProtonLike( const caf::SRSliceProxy* slice, const unsigned int idx ) {
    auto const& trk = slice->reco.pfp.at(idx).trk; //GetTrack( slc, idxTrk );
    if ( std::isnan(slice->vertex.x) || std::isnan(slice->vertex.y) || std::isnan(slice->vertex.z) ) return false;
    const float Atslc = std::hypot(slice->vertex.x - trk.start.x,
                                  slice->vertex.y - trk.start.y,
                                  slice->vertex.z - trk.start.z);
    const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slice,idx));
    if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) return false;
    const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
    const float Chi2Proton = trk.chi2pid[2].chi2_proton;
    const float Chi2Muon = trk.chi2pid[2].chi2_muon;
    if ( !(Contained && Chi2Proton < 90. && Chi2Muon > 30. ) ) return false;
    return true;
  }

  bool IsMuonLike( const caf::SRSliceProxy* slice, const unsigned int idx ) {
    auto const& trk = slice->reco.pfp.at(idx).trk; //GetTrack( slc, idxTrk );
    
    if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) return false;
    if ( std::isnan(slice->vertex.x) || std::isnan(slice->vertex.y) || std::isnan(slice->vertex.z) ) return false;
    const float Atslc = std::hypot( slice->vertex.x - trk.start.x,
                                    slice->vertex.y - trk.start.y,
                                    slice->vertex.z - trk.start.z);
    const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slice,idx));

    if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) return false;
    const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
    
    const float Chi2Proton = trk.chi2pid[2].chi2_proton;
    const float Chi2Muon = trk.chi2pid[2].chi2_muon;
    //track len of muon
    if ( (!Contained && trk.len > 5.) || ( Contained && trk.len > 5. && Chi2Proton > 60. && Chi2Muon < 30.) ) return true;
    return false;
  }

  bool IsPionLike( const caf::SRSliceProxy* slice, const unsigned int idx ) {
    auto const& trk = slice->reco.pfp.at(idx).trk; //GetTrack( slc, idxTrk );
    if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) return false;
    if ( std::isnan(slice->vertex.x) || std::isnan(slice->vertex.y) || std::isnan(slice->vertex.z) ) return false;
    const float Atslc = std::hypot( slice->vertex.x - trk.start.x,
                                    slice->vertex.y - trk.start.y,
                                    slice->vertex.z - trk.start.z);
    const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slice,idx));

    if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) return false;
    const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
    
    const float Chi2Proton = trk.chi2pid[2].chi2_proton;
    const float Chi2Muon = trk.chi2pid[2].chi2_muon;
    //track len of muon
    if ( Contained && Chi2Proton > 60. && Chi2Muon < 30.) return true;
    return false;
  }

  bool IsPrimaryPFP( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
    return slice->reco.pfp.at(idxTrk).parent_is_primary;
  }

  // Dummy vars
  const Var kNuMIDummyVar1([](const caf::SRSliceProxy* slc) -> int {
    return 1;
  });

  const Var kNuMIDummyVar0([](const caf::SRSliceProxy* slc) -> int {
    return 0;
  });

  // Emulated trigger time
  const SpillVar kNuMISpillTriggerTime ( [](const caf::SRSpillProxy *sr) -> double {
    double triggerTime = 0.;

    double foundTriggerTime = sr->hdr.triggerinfo.trigger_within_gate;
    if ( !std::isnan(foundTriggerTime) && !std::isinf(foundTriggerTime) && foundTriggerTime < -15. ) triggerTime = -15.;
    else if ( std::isnan(foundTriggerTime) || std::isinf(foundTriggerTime) ) triggerTime = -16.;
    else if ( foundTriggerTime > 30. ) triggerTime = 30.; //1.6 for BNB
    else triggerTime = foundTriggerTime;

    return triggerTime;
  });

  //const SpillVar kCRTMPMT_FlashMatching( [](const caf::SRSpillProxy* sr) -> float {
  //  const auto& match = sr->crtpmt_matches;
  //  double flashGateTime = match.flashGateTime;
  //  return flashGateTime;
  //});

  // Muon candidate
  const Var kNuMIMuonCandidateIdx([](const caf::SRSliceProxy* slc) -> int {
      float Longest(0);
      int PTrackInd(-3);
      float p(-5.f);

      unsigned int idxTrk = 0;
      while ( IsValidTrkIdx(slc, idxTrk) ) {
        if ( !IsTracklikeTrack(slc, idxTrk) ) { 
          idxTrk+=1;
          continue;
        }
        auto const& trk = slc->reco.pfp.at(idxTrk).trk;
        unsigned int thisIdx = idxTrk;
        idxTrk+=1;

        if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
        if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return -1;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

        if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) continue;
        const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);

        if(Contained) {
          p = trk.rangeP.p_muon;
        }else{
          if(std::isnan(trk.mcsP.fwdP_muon) || std::isinf(trk.mcsP.fwdP_muon) )  p = -5.0;
          else p = trk.mcsP.fwdP_muon;
        }
        
        if ( p < 0.226 ) continue;
        
        const float Chi2Proton = trk.chi2pid[2].chi2_proton;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        //track len of muon
        if (     (!Contained && trk.len > 5.) 
              || ( Contained && trk.len > 5. && Chi2Proton > 60. && Chi2Muon < 30.) ) {
          if ( trk.len <= Longest ) continue;
          Longest = trk.len;
          PTrackInd = thisIdx;
        }
      }

      return PTrackInd;
    });

  // Proton candidate
  const Var kNuMIProtonCandidateIdx([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kNuMIMuonCandidateIdx(slc);

    float Longest(0);
    int PTrackInd(-2);

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      if ( !IsTracklikeTrack(slc, idxTrk) ) { 
        idxTrk+=1;
        continue;
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk; //GetTrack( slc, idxTrk );
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return -1;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

      if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) continue;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( Contained && Chi2Proton < 90. && Chi2Muon > 30. ) {
        if ( trk.len <= Longest ) continue;
        Longest = trk.len;
        PTrackInd = thisIdx;
      }
    }

    return PTrackInd;
  });

  // MultiVar for the charged pion candidate indices
  const MultiVar kNuMIChargedPionCandidateIdxs([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;

    int primaryInd = kNuMIMuonCandidateIdx(slc);
    int primaryProtonInd = kNuMIProtonCandidateIdx(slc);

    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); ++i_pfp){

      if ( i_pfp == (unsigned int)primaryInd || i_pfp == (unsigned int)primaryProtonInd ) {
        continue; // skip the particle which is the muon or leading proton candidate!
      }

      auto const& trk = slc->reco.pfp.at(i_pfp).trk;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) continue;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,i_pfp));

      if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) continue;

      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( Contained && Chi2Proton > 60. && Chi2Muon < 30. ){
        rets.push_back( i_pfp );
      }
    }

    return rets;

  });

  // MultiVar for the photon candidate indices from reconstruction
  const MultiVar kNuMIPhotonCandidateIdxs([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;

    int primaryInd = kNuMIMuonCandidateIdx(slc);
    //int primaryProtonInd = kNuMIProtonCandidateIdx(slc);

    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); ++i_pfp){

      if (     i_pfp == (unsigned int)primaryInd 
            //|| i_pfp == (unsigned int)primaryProtonInd 
          ) continue; // skip the particle which is the muon or leading proton candidate!
      
      //if ( !IsShowerlike(slc, i_pfp) ) { 
      //  continue; // skip things with track score < 0.45
      //}
      auto const& shw = slc->reco.pfp.at(i_pfp).shw;
      //auto const& trk = slc->reco.pfp.at(i_pfp).trk;

      // Check if shower fit even seems kind-of valid:
      if (    std::isnan(shw.start.x)
           //|| (shw.start.x > -5.5 && shw.start.x < -4.5)
           || std::isnan(shw.len) || shw.len <= 0. ) continue;
      
      // if it meets this then we're not going to cut on it...
      if ( std::isnan(slc->reco.pfp.at(i_pfp).trackScore) || std::isinf(slc->reco.pfp.at(i_pfp).trackScore) || slc->reco.pfp.at(i_pfp).trackScore <= 0. ) continue;
      //if ( std::isnan(shw.plane[2].energy) || std::isinf(shw.plane[2].energy) || shw.plane[2].energy <= 0.02 ) continue;
      if ( std::isnan(shw.bestplane_energy) || std::isinf(shw.bestplane_energy) || shw.bestplane_energy <= 0.00 ) continue;

      // and... if it meets then then we're not going to cut on it...
      if ( std::isnan(shw.conversion_gap) || std::isinf(shw.conversion_gap) || shw.conversion_gap <= 0. ) continue;
      if ( !isInFV(shw.start.x,shw.start.y,shw.start.z) ) continue; //Quality cut
      //if ( IsProtonLike(slc,i_pfp) ) continue; //Smells like a proton
      //if ( IsMuonLike(slc,i_pfp) ) continue; //Smells like a muon
      //if ( IsPionLike(slc,i_pfp) ) continue; //Smells like a pion
      
      //if ( std::isnan(trk.chi2pid[2].chi2_proton) || std::isinf(trk.chi2pid[2].chi2_proton) || trk.chi2pid[2].chi2_proton < 120. ) continue; //Cut to reduce protons
      // if we got here, then it should be the case that the fit seems valid and:
      // shwE > 0.040 GeV
      // trackScore < 0.45 (technically <= 0.45)
      // conversionGap > 5. cm

      rets.push_back( i_pfp );

    }

    // guess we're not cutting anything
    return rets;

  });

  // Number of truth-matched photons in reconstruction
  const Var kNumberTruthMatchRecoPhotons([](const caf::SRSliceProxy* slc) -> int {
    int nphotons(0);

    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); ++i_pfp){
      auto const& shw = slc->reco.pfp.at(i_pfp).shw;
      
      // Get truth matched shower pdg
      if ( std::isnan(shw.truth.p.pdg) || std::isinf(shw.truth.p.pdg) ) continue;
      if ( abs(shw.truth.p.pdg) != 22 ) continue;
      if ( shw.truth.p.start_process != 3 ) continue;
      nphotons++;
    }
    return nphotons;
  });

  const Var kNumberTruthMatchRecoMuons([](const caf::SRSliceProxy* slc) -> int {
    int nmuons(0);
    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); ++i_pfp){
      auto const& trk = slc->reco.pfp.at(i_pfp).trk;
      
      // Get truth matched shower pdg
      if ( std::isnan(trk.truth.p.pdg) || std::isinf(trk.truth.p.pdg) ) continue;
      if ( abs(trk.truth.p.pdg) != 13 ) continue;
      if ( trk.truth.p.start_process != 0 ) continue;
      nmuons++;
    }
    return nmuons;
  });

  // Leading photon candidate
  const Var kNuMILeadingPhotonCandidateIdx([](const caf::SRSliceProxy* slc) -> int {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if ( photon_indices.size() == 0 ) return -5;

    // find the leading photon candidate
    int leadingPhotonIdx = -1;
    double leadingPhotonE = -1.;

    for ( auto const& idx : photon_indices ) {
      auto const& shw = slc->reco.pfp.at(idx).shw;
      
      //auto const& trk = slc->reco.pfp.at(idx).trk;
      //Quality Cuts
      //if ( std::isnan(trk.chi2pid[2].chi2_proton) || std::isinf(trk.chi2pid[2].chi2_proton) || trk.chi2pid[2].chi2_proton < 90. ) continue; //Cut to reduce protons
      //if ( shw.plane[2].energy < 0.075 ) continue; //Quality cut
      //if ( std::isnan(slc->reco.pfp[idx].trackScore) || std::isinf(slc->reco.pfp[idx].trackScore) || slc->reco.pfp[idx].trackScore > 0.55 ) continue; //Quality cut trackScore float trkscore = slc->reco.pfp[idx].trackScore;
      //if ( shw.plane[2].nHits < 20 ) continue; //Quality cut
     

      // Find shower with highest energy
      if ( shw.bestplane_energy > leadingPhotonE ) {
        //leadingPhotonE = shw.plane[2].energy;
        leadingPhotonE = shw.bestplane_energy;
        leadingPhotonIdx = idx;
      }
    }
    return leadingPhotonIdx;
  });

  // SubLeading photon candidate
  const Var kNuMISubLeadingPhotonCandidateIdx([](const caf::SRSliceProxy* slc) -> int {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if ( photon_indices.size() <= 1 ) return -5;
    int leadingPhotonIdx = kNuMILeadingPhotonCandidateIdx(slc);
    
    // find the subleading photon candidate
    int subleadingPhotonIdx = -5;
    double subleadingPhotonE = -5.;
    for ( auto const& idx : photon_indices ) {
      if ( idx == leadingPhotonIdx ) continue;
      auto const& shw = slc->reco.pfp.at(idx).shw;
      //auto const& trk = slc->reco.pfp.at(idx).trk;

      //Quality Cuts
      //if ( shw.plane[2].energy < 0.02 ) continue; //Quality cut
      //if ( std::isnan(slc->reco.pfp[idx].trackScore) || std::isinf(slc->reco.pfp[idx].trackScore) || slc->reco.pfp[idx].trackScore > 0.6 ) continue; //Quality cut trackScore float trkscore = slc->reco.pfp[idx].trackScore;
      
      // Find shower with second highest energy
      if ( shw.bestplane_energy > subleadingPhotonE ) {
        //subleadingPhotonE = shw.plane[2].energy;
        subleadingPhotonE = shw.bestplane_energy;
        subleadingPhotonIdx = idx;
      }
    }
    return subleadingPhotonIdx;
  });

  const Var kNumberRecoShowers([](const caf::SRSliceProxy* slc) -> int {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc); //track score > 0.45
    if (std::isnan(photon_indices.size())) return 0;
    int numshowers = photon_indices.size();
    return numshowers;
  });

  const Var kNumberPFPs([](const caf::SRSliceProxy* slc) -> int {
    return slc->reco.npfp;
  });

  const Var kNumberChargedPions([](const caf::SRSliceProxy* slc) -> int {
    std::vector<double> chargedpion_indices = kNuMIChargedPionCandidateIdxs(slc);
    if (std::isnan(chargedpion_indices.size())) return 0;
    int numchargedpions = chargedpion_indices.size();
    return numchargedpions;
  });

  // MultiVar for the track candidate indices
  const MultiVar kNuMITrackCandidateIdxs([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;

    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); ++i_pfp){

      if ( !IsTracklikeTrack(slc, i_pfp) ) { 
        continue; // skip things with track score <= 0.45
      }
      auto const& trk = slc->reco.pfp.at(i_pfp).trk;

      // Check if track fit even seems kind-of valid:
      if ( std::isnan(trk.start.x) || (trk.start.x > -5.5 && trk.start.x < -4.5) ||
           std::isnan(trk.len) || trk.len <= 0. ) continue;

      rets.push_back( i_pfp );

    }

    // guess we're not cutting anything
    return rets;

  });

  const Var kNumberRecoTracks([](const caf::SRSliceProxy* slc) -> int {
    std::vector<double> trk_indices = kNuMITrackCandidateIdxs(slc); //track score <= 0.45
    if (std::isnan(trk_indices.size())) return 0;
    int numtracks = trk_indices.size();
    return numtracks;
  });


  // Reco muon momentum
  const Var kNuMIMuonCandidateRecoP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kNuMIMuonCandidateIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      if(std::isnan(trk.mcsP.fwdP_muon)) return p = -5.0;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained) p = trk.rangeP.p_muon;
      //else p = trk.mcsP.fwdP_muon;
      else p = 1;
    }
    return p;
  });

  // True muon momentum
  const Var kNuMIMuonTrueP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    if ( slc->truth.index >= 0 ) p = PrimaryUtil::MuonP_True(slc->truth);

    return p;
  });

  // Reco muon track length
  const Var kNuMIRecoMuonLength([](const caf::SRSliceProxy* slc) -> double {
    float ret(-5.f);
    int candIdx = kNuMIMuonCandidateIdx(slc);
    if( candIdx >= 0 ){
      auto const& trk = slc->reco.pfp.at(candIdx).trk;

      ret = trk.len;
    }

     return ret;
  });
  // True muon length
  const Var kNuMITrueMuonLength([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    if ( slc->truth.index >= 0 ) p = PrimaryUtil::MuonLength_True(slc->truth);

    return p;
  });

  // Reco muon transverse momentum
  const Var kNuMIRecoMuonPt([](const caf::SRSliceProxy* slc) -> double {
    float ret(-5.f);
    int candIdx = kNuMIMuonCandidateIdx(slc);
    if( candIdx >= 0 ){
      auto const& trk = slc->reco.pfp.at(candIdx).trk;
      double momentum = kNuMIMuonCandidateRecoP(slc);
      TVector3 vec_momentum(trk.dir.x, trk.dir.y, trk.dir.z);
      vec_momentum *= momentum;

      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      dFromNuMI *= 100.; // m to cm
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

      TVector3 vec_pt = vec_momentum - (vec_momentum.Dot(vec_numi_to_vtx)) * vec_numi_to_vtx;

      ret = vec_pt.Mag();

    }

     return ret;
  });
  // True muon transverse momentum
  const Var kNuMITrueMuonPt([](const caf::SRSliceProxy* slc) -> double {
    float ret(-5.f);
    if ( slc->truth.index >= 0 ) ret = PrimaryUtil::MuonPt_True(slc->truth);

    return ret;
  });

  // Reco proton momentum
  const Var kNuMIProtonCandidateRecoP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    if ( kNuMIProtonCandidateIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kNuMIProtonCandidateIdx(slc)).trk;
      p = trk.rangeP.p_proton;
    }
    return p;
  });

  // True proton momentum
  const Var kNuMIProtonTrueP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    if ( slc->truth.index >= 0 ) p = PrimaryUtil::ProtonP_True(slc->truth);

    return p;
  });

  // Reco proton transverse momentum
  const Var kNuMIRecoProtonPt([](const caf::SRSliceProxy* slc) -> double {
    float ret(-5.f);
    int candIdx = kNuMIProtonCandidateIdx(slc);
    if( candIdx >= 0 ){
      auto const& trk = slc->reco.pfp.at(candIdx).trk;
      double momentum = kNuMIProtonCandidateRecoP(slc);
      TVector3 vec_momentum(trk.dir.x, trk.dir.y, trk.dir.z);
      vec_momentum *= momentum;

      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      dFromNuMI *= 100.; // m to cm
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

      TVector3 vec_pt = vec_momentum - (vec_momentum.Dot(vec_numi_to_vtx)) * vec_numi_to_vtx;

      ret = vec_pt.Mag();

    }

     return ret;
  });
  // True proton transverse momentum
  const Var kNuMITrueProtonPt([](const caf::SRSliceProxy* slc) -> double {
    float ret(-5.f);
    if ( slc->truth.index >= 0 ) ret = PrimaryUtil::ProtonPt_True(slc->truth);

    return ret;
  });

  // Reco proton track length
  const Var kNuMIRecoProtonLength([](const caf::SRSliceProxy* slc) -> double {
    float ret(-5.f);
    int candIdx = kNuMIProtonCandidateIdx(slc);
    if( candIdx >= 0 ){
      auto const& trk = slc->reco.pfp.at(candIdx).trk;

      ret = trk.len;
    }

     return ret;
  });
  // True proton length
  const Var kNuMITrueProtonLength([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    if ( slc->truth.index >= 0 ) p = PrimaryUtil::ProtonLength_True(slc->truth);

    return p;
  });

  // Reco CosTh(numi)
  const Var kNuMIRecoCosThBeam([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if ( kNuMIMuonCandidateIdx(slc) >= 0 ) {
      auto const& mutrk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;

      double magNuMI = sqrt(315.120380*315.120380 + 33.644912*33.644912 + 733.632532*733.632532);
      TVector3 rFromNuMI(315.120380/magNuMI, 33.644912/magNuMI, 733.632532/magNuMI);

      TVector3 muDir(mutrk.dir.x, mutrk.dir.y, mutrk.dir.z);
      muDir = muDir.Unit();

      costh = TMath::Cos( muDir.Angle(rFromNuMI) );
    }

    return costh;
  });

  // True CosTh(numi)
  const Var kNuMITrueCosThBeam([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    if ( slc->truth.index >= 0 ) costh = PrimaryUtil::MuonCosThBeam_True(slc->truth);

    return costh;
  });

  // Truth angle between photons from neutral pion decay
  const Var kNuMITrueCosThPhotonPhoton([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    if ( slc->truth.index >= 0 ) costh = PrimaryUtil::CosThPhotonPhoton_True(slc->truth);
    return costh;
  });

  // Leading true photon G4ID
  const Var kNuMITrueLeadingPhotonG4ID([](const caf::SRSliceProxy* slc) -> int {
    int g4id(-5);
    if ( slc->truth.index >= 0 ) g4id = PrimaryUtil::Pi0LeadingPhotonG4ID(slc->truth);
    return g4id;
  });

  // Subleading true photon G4ID
  const Var kNuMITrueSubLeadingPhotonG4ID([](const caf::SRSliceProxy* slc) -> int {
    int g4id(-5);
    if ( slc->truth.index >= 0 ) g4id = PrimaryUtil::Pi0SubLeadingPhotonG4ID(slc->truth);
    return g4id;
  });

  // Reco Muon angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  const Var kNuMIRecoCosThVtx([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    int candIdx = kNuMIMuonCandidateIdx(slc);
    if ( candIdx>=0 ) {
      auto const& trk = slc->reco.pfp.at(candIdx).trk;
      TVector3 vec_trk(trk.dir.x, trk.dir.y, trk.dir.z);

      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      dFromNuMI *= 100.; // m to cm
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

      double angle = vec_trk.Angle(vec_numi_to_vtx);
      costh = TMath::Cos(angle);
    }

    return costh;
  });
  // Reco Muon angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  const Var kNuMITrueCosThVtx([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    if ( slc->truth.index >= 0 ) costh = PrimaryUtil::MuonNuCosineTheta_True(slc->truth);

    return costh;
  });

  // Reco Proton angle w.r.t. beam
  const Var kNuMIProtonRecoCosThBeam([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    int candIdx = kNuMIProtonCandidateIdx(slc);
    if ( candIdx >= 0 ) {
      auto const& trk = slc->reco.pfp.at(candIdx).trk;

      double magNuMI = sqrt(315.120380*315.120380 + 33.644912*33.644912 + 733.632532*733.632532);
      TVector3 rFromNuMI(315.120380/magNuMI, 33.644912/magNuMI, 733.632532/magNuMI);

      TVector3 trkDir(trk.dir.x, trk.dir.y, trk.dir.z);
      trkDir = trkDir.Unit();

      costh = TMath::Cos( trkDir.Angle(rFromNuMI) );
    }

    return costh;
  });
  // True Proton angle w.r.t. beam
  const Var kNuMIProtonTrueCosThBeam([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    if ( slc->truth.index >= 0 ) costh = PrimaryUtil::ProtonCosThBeam_True(slc->truth);

    return costh;
  });

  // Reco Proton angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  const Var kNuMIProtonRecoCosThVtx([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    int candIdx = kNuMIProtonCandidateIdx(slc);
    if ( candIdx>=0 ) {
      auto const& trk = slc->reco.pfp.at(candIdx).trk;
      TVector3 vec_trk(trk.dir.x, trk.dir.y, trk.dir.z);

      const auto& vtx = slc->vertex;
      TVector3 vec_vtx(vtx.x, vtx.y, vtx.z);
      TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      dFromNuMI *= 100.; // m to cm
      TVector3 vec_numi_to_vtx = (dFromNuMI + vec_vtx).Unit();

      double angle = vec_trk.Angle(vec_numi_to_vtx);
      costh = TMath::Cos(angle);
    }

    return costh;
  });
  // Reco Proton angle w.r.t. numi-to-vtx direction (= proxy of neutrino direction)
  const Var kNuMIProtonTrueCosThVtx([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    if ( slc->truth.index >= 0 ) costh = PrimaryUtil::ProtonNuCosineTheta_True(slc->truth);

    return costh;
  });


  // Reco CosTh(mu,p)
  const Var kNuMIRecoCosThMuP([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if ( kNuMIMuonCandidateIdx(slc) >= 0 && kNuMIProtonCandidateIdx(slc) >= 0 ) {
      auto const& mutrk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      auto const& ptrk = slc->reco.pfp.at(kNuMIProtonCandidateIdx(slc)).trk;

      TVector3 muDir(mutrk.dir.x, mutrk.dir.y, mutrk.dir.z);
      muDir = muDir.Unit();
      TVector3 pDir(ptrk.dir.x, ptrk.dir.y, ptrk.dir.z);
      pDir = pDir.Unit();

      costh = TMath::Cos( muDir.Angle(pDir) );
    }

    return costh;
  });

  // True CosTh(mu,p)
  const Var kNuMITrueCosThMuP([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);
    if ( slc->truth.index >= 0 ) costh = PrimaryUtil::CosThMuonProton_True(slc->truth);

    return costh;
  });

  const Var kNuMIRecodeltaPT([](const caf::SRSliceProxy* slc) -> float {
    float ret(-5.f);

    int MuonTrackIndex = kNuMIMuonCandidateIdx(slc);
    int ProtonTrackIndex = kNuMIProtonCandidateIdx(slc);
    if ( MuonTrackIndex >= 0 && ProtonTrackIndex >= 0 ) {

      auto const& mu_trk = slc->reco.pfp.at(MuonTrackIndex).trk;
      double mu_trk_p = kNuMIMuonCandidateRecoP(slc);
      TVector3 vec_p_mu(mu_trk.dir.x, mu_trk.dir.y, mu_trk.dir.z);
      vec_p_mu *= mu_trk_p;

      auto const& pro_trk = slc->reco.pfp.at(ProtonTrackIndex).trk;
      double pro_trk_p = kNuMIProtonCandidateRecoP(slc);
      TVector3 vec_p_pro(pro_trk.dir.x, pro_trk.dir.y, pro_trk.dir.z);
      vec_p_pro *= pro_trk_p;

      const auto& vtx = slc->vertex;
      // ICARUSZero-to-vertex
      TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
      // NuMI-to-ICARUSZero
      TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);
      vec_NuMI_to_ICARUS *= 100.; // meter to centimeter
      // NuMI-to-vertex = NuMI-to-ICARUSZero + ICARUSZero-to-vertex
      TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

      ret = PrimaryUtil::CalcTKI_deltaPT(vec_p_mu, vec_p_pro, unit_numi_to_vtx);

    }

    return ret;

  });

  const Var kNuMITruedeltaPT([](const caf::SRSliceProxy* slc) -> float {
    float ret(-5.f);
    if ( slc->truth.index >= 0 ) ret = PrimaryUtil::deltaPT_True(slc->truth);

    return ret;
  });

  const Var kNuMIRecodeltaPTx([](const caf::SRSliceProxy* slc) -> float {
    float ret(-99999); // TODO deltaPTx is a signed variable.. for now just making it very very largly negative

    int MuonTrackIndex = kNuMIMuonCandidateIdx(slc);
    int ProtonTrackIndex = kNuMIProtonCandidateIdx(slc);
    if ( MuonTrackIndex >= 0 && ProtonTrackIndex >= 0 ) {

      auto const& mu_trk = slc->reco.pfp.at(MuonTrackIndex).trk;
      double mu_trk_p = kNuMIMuonCandidateRecoP(slc);
      TVector3 vec_p_mu(mu_trk.dir.x, mu_trk.dir.y, mu_trk.dir.z);
      vec_p_mu *= mu_trk_p;

      auto const& pro_trk = slc->reco.pfp.at(ProtonTrackIndex).trk;
      double pro_trk_p = kNuMIProtonCandidateRecoP(slc);
      TVector3 vec_p_pro(pro_trk.dir.x, pro_trk.dir.y, pro_trk.dir.z);
      vec_p_pro *= pro_trk_p;

      const auto& vtx = slc->vertex;
      // ICARUSZero-to-vertex
      TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
      // NuMI-to-ICARUSZero
      TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);
      vec_NuMI_to_ICARUS *= 100.; // meter to centimeter
      // NuMI-to-vertex = NuMI-to-ICARUSZero + ICARUSZero-to-vertex
      TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

      ret =  PrimaryUtil::CalcTKI_deltaPTx(vec_p_mu, vec_p_pro, unit_numi_to_vtx);

    }

    return ret;

  });

  const Var kNuMITruedeltaPTx([](const caf::SRSliceProxy* slc) -> float {
    float ret(-99999); // TODO deltaPTx is a signed variable.. for now just making it very very largly negative
    if ( slc->truth.index >= 0 ) ret = PrimaryUtil::deltaPTx_True(slc->truth);

    return ret;
  });

  const Var kNuMIRecodeltaPTy([](const caf::SRSliceProxy* slc) -> float {
    float ret(-99999); // TODO deltaPTy is a signed variable.. for now just making it very very largly negative

    int MuonTrackIndex = kNuMIMuonCandidateIdx(slc);
    int ProtonTrackIndex = kNuMIProtonCandidateIdx(slc);
    if ( MuonTrackIndex >= 0 && ProtonTrackIndex >= 0 ) {

      auto const& mu_trk = slc->reco.pfp.at(MuonTrackIndex).trk;
      double mu_trk_p = kNuMIMuonCandidateRecoP(slc);
      TVector3 vec_p_mu(mu_trk.dir.x, mu_trk.dir.y, mu_trk.dir.z);
      vec_p_mu *= mu_trk_p;

      auto const& pro_trk = slc->reco.pfp.at(ProtonTrackIndex).trk;
      double pro_trk_p = kNuMIProtonCandidateRecoP(slc);
      TVector3 vec_p_pro(pro_trk.dir.x, pro_trk.dir.y, pro_trk.dir.z);
      vec_p_pro *= pro_trk_p;

      const auto& vtx = slc->vertex;
      // ICARUSZero-to-vertex
      TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
      // NuMI-to-ICARUSZero
      TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);
      vec_NuMI_to_ICARUS *= 100.; // meter to centimeter
      // NuMI-to-vertex = NuMI-to-ICARUSZero + ICARUSZero-to-vertex
      TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

      ret = PrimaryUtil::CalcTKI_deltaPTy(vec_p_mu, vec_p_pro, unit_numi_to_vtx);

    }

    return ret;
  });

  const Var kNuMITruedeltaPTy([](const caf::SRSliceProxy* slc) -> float {
    float ret(-99999); // TODO deltaPTy is a signed variable.. for now just making it very very largly negative
    if ( slc->truth.index >= 0 ) ret = PrimaryUtil::deltaPTy_True(slc->truth);

    return ret;
  });

  const Var kNuMIRecodeltaalphaT([](const caf::SRSliceProxy* slc) -> float {
    float ret(-5.f); // acos is in the interval of [0,pi]

    int MuonTrackIndex = kNuMIMuonCandidateIdx(slc);
    int ProtonTrackIndex = kNuMIProtonCandidateIdx(slc);
    if ( MuonTrackIndex >= 0 && ProtonTrackIndex >= 0 ) {

      auto const& mu_trk = slc->reco.pfp.at(MuonTrackIndex).trk;
      double mu_trk_p = kNuMIMuonCandidateRecoP(slc);
      TVector3 vec_p_mu(mu_trk.dir.x, mu_trk.dir.y, mu_trk.dir.z);
      vec_p_mu *= mu_trk_p;

      auto const& pro_trk = slc->reco.pfp.at(ProtonTrackIndex).trk;
      double pro_trk_p = kNuMIProtonCandidateRecoP(slc);
      TVector3 vec_p_pro(pro_trk.dir.x, pro_trk.dir.y, pro_trk.dir.z);
      vec_p_pro *= pro_trk_p;

      const auto& vtx = slc->vertex;
      // ICARUSZero-to-vertex
      TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
      // NuMI-to-ICARUSZero
      TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);
      vec_NuMI_to_ICARUS *= 100.; // meter to centimeter
      // NuMI-to-vertex = NuMI-to-ICARUSZero + ICARUSZero-to-vertex
      TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

      ret = PrimaryUtil::CalcTKI_deltaalphaT(vec_p_mu, vec_p_pro, unit_numi_to_vtx);

    }

    return ret;
  });

  const Var kNuMITruedeltaalphaT([](const caf::SRSliceProxy* slc) -> float {
    float ret(-5.f); // acos is in the interval of [0,pi]
    if ( slc->truth.index >= 0 ) ret = PrimaryUtil::deltaalphaT_True(slc->truth);

    return ret;
  });

  const Var kNuMIRecodeltaphiT([](const caf::SRSliceProxy* slc) -> float {
    float ret(-5.f); // acos is in the interval of [0,pi]

    int MuonTrackIndex = kNuMIMuonCandidateIdx(slc);
    int ProtonTrackIndex = kNuMIProtonCandidateIdx(slc);
    if ( MuonTrackIndex >= 0 && ProtonTrackIndex >= 0 ) {

      auto const& mu_trk = slc->reco.pfp.at(MuonTrackIndex).trk;
      double mu_trk_p = kNuMIMuonCandidateRecoP(slc);
      TVector3 vec_p_mu(mu_trk.dir.x, mu_trk.dir.y, mu_trk.dir.z);
      vec_p_mu *= mu_trk_p;

      auto const& pro_trk = slc->reco.pfp.at(ProtonTrackIndex).trk;
      double pro_trk_p = kNuMIProtonCandidateRecoP(slc);
      TVector3 vec_p_pro(pro_trk.dir.x, pro_trk.dir.y, pro_trk.dir.z);
      vec_p_pro *= pro_trk_p;

      const auto& vtx = slc->vertex;
      // ICARUSZero-to-vertex
      TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
      // NuMI-to-ICARUSZero
      TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);
      vec_NuMI_to_ICARUS *= 100.; // meter to centimeter
      // NuMI-to-vertex = NuMI-to-ICARUSZero + ICARUSZero-to-vertex
      TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

      ret = PrimaryUtil::CalcTKI_deltaphiT(vec_p_mu, vec_p_pro, unit_numi_to_vtx);

    }

    return ret;
  });

  const Var kNuMITruedeltaphiT([](const caf::SRSliceProxy* slc) -> float {
    float ret(-5.f); // acos is in the interval of [0,pi]
    if ( slc->truth.index >= 0 ) ret = PrimaryUtil::deltaphiT_True(slc->truth);

    return ret;
  });



  //////////////////////////////
  //Reco Nu interaction Calcualtions
  //////////////////////////////

  const Var kSlcVertexX([](const caf::SRSliceProxy* slc) -> float {
    if (std::isnan(slc->vertex.x) || std::isinf(slc->vertex.x)) {
      return -9999.f;
    }
    return slc->vertex.x;
  });

  const Var kSlcVertexY([](const caf::SRSliceProxy* slc) -> float {
    if (std::isnan(slc->vertex.y) || std::isinf(slc->vertex.y)) {
      return -9999.f;
    }
    return slc->vertex.y;
  });

  const Var kSlcVertexZ([](const caf::SRSliceProxy* slc) -> float {
    if (std::isnan(slc->vertex.z) || std::isinf(slc->vertex.z)) {
      return -9999.f;
    }
    return slc->vertex.z;
  });

  //////////////////////////////
  //Reco Muon Calcualtions
  //////////////////////////////
  const Var KMuonCandidateRecoStartX([](const caf::SRSliceProxy* slc) -> float {
    float ret(-9999.f);
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return ret;
    auto const& trk = slc->reco.pfp.at(Idx).trk;
    ret = trk.start.x;
    return ret;
  });

  const Var KMuonCandidateRecoStartY([](const caf::SRSliceProxy* slc) -> float {
    float ret(-9999.f);
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return ret;
    auto const& trk = slc->reco.pfp.at(Idx).trk;
    ret = trk.start.y;
    return ret;
  });

  const Var KMuonCandidateRecoStartZ([](const caf::SRSliceProxy* slc) -> float {
    float ret(-9999.f);
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return ret;
    auto const& trk = slc->reco.pfp.at(Idx).trk;
    ret = trk.start.z;
    return ret;
  });

  const Var kMuonCandidateTrueStartX([](const caf::SRSliceProxy* slc) -> float {
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return -9999.f;
    if ( std::isnan(slc->reco.pfp.at(Idx).trk.truth.p.start.x) || std::isinf(slc->reco.pfp.at(Idx).trk.truth.p.start.x) ) return -9999.f;
    return slc->reco.pfp.at(Idx).trk.truth.p.start.x;
  });

  const Var kMuonCandidateTrueStartY([](const caf::SRSliceProxy* slc) -> float {
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return -9999.f;
    if ( std::isnan(slc->reco.pfp.at(Idx).trk.truth.p.start.y) || std::isinf(slc->reco.pfp.at(Idx).trk.truth.p.start.y) ) return -9999.f;
    return slc->reco.pfp.at(Idx).trk.truth.p.start.y;
  });

  const Var kMuonCandidateTrueStartZ([](const caf::SRSliceProxy* slc) -> float {
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return -9999.f;
    if ( std::isnan(slc->reco.pfp.at(Idx).trk.truth.p.start.z) || std::isinf(slc->reco.pfp.at(Idx).trk.truth.p.start.z) ) return -9999.f;
    return slc->reco.pfp.at(Idx).trk.truth.p.start.z;
  });
  
  const Var kMuonCandidatePDG([](const caf::SRSliceProxy* slc) -> float {
    int Idx = kNuMIMuonCandidateIdx(slc);
    if (Idx < 0) return -9999.f;
    if ( std::isnan(slc->reco.pfp.at(Idx).trk.truth.p.pdg) || std::isinf(slc->reco.pfp.at(Idx).trk.truth.p.pdg) ) return -9999.f;
    return slc->reco.pfp.at(Idx).trk.truth.p.pdg;
  });

  const Var kProtonCandidatePDG([](const caf::SRSliceProxy* slc) -> float {
    int Idx = kNuMIProtonCandidateIdx(slc);
    if (Idx < 0) return -9999.f;
    if ( std::isnan(slc->reco.pfp.at(Idx).trk.truth.p.pdg) || std::isinf(slc->reco.pfp.at(Idx).trk.truth.p.pdg) ) return -9999.f;
    return slc->reco.pfp.at(Idx).trk.truth.p.pdg;
  });

  /////////////////////////////////////////////
  /////Photon Calculations start here//////////
  /////////////////////////////////////////////

  const Var kNuMILeadingPhotonCandidateE([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    //return slc->reco.pfp[idx].shw.plane[2].energy;//shw.plane[2].energy
    return slc->reco.pfp[idx].shw.bestplane_energy;
  });

  const Var kNuMILeadingPhotonCandidateTrueE([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.genE) || std::isinf(slc->reco.pfp[idx].shw.truth.p.genE) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.genE;
  });

  const Var kNuMISubLeadingPhotonCandidateE([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    //return slc->reco.pfp[idx].shw.plane[2].energy;
    return slc->reco.pfp[idx].shw.bestplane_energy;
  });

  const Var kNuMISubLeadingPhotonCandidateTrueE([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.genE) || std::isinf(slc->reco.pfp[idx].shw.truth.p.genE) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.genE;
  });

  const Var kNuMIPhotonCandidatesOpeningAngle([](const caf::SRSliceProxy* slc) -> float {
    int idxMaxE = kNuMILeadingPhotonCandidateIdx(slc);
    int idxScdy = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idxMaxE<0 || idxScdy<0) return -5.f;

    // Now get the quantities we want from this, following results of Jamie's w.r.t. the opening angle:
    // Uses vertex to shower:
    auto const& shw1 = slc->reco.pfp[idxMaxE].shw;
    TVector3 vecShw1(shw1.start.x - slc->vertex.x, shw1.start.y - slc->vertex.y, shw1.start.z - slc->vertex.z);
    auto const& shw2 = slc->reco.pfp[idxScdy].shw;
    TVector3 vecShw2(shw2.start.x - slc->vertex.x, shw2.start.y - slc->vertex.y, shw2.start.z - slc->vertex.z);
    float openAngle = vecShw1.Angle(vecShw2);

    if (std::isinf(openAngle) || std::isnan(openAngle)) return -5.f;

    return openAngle;
  });

  const Var kNuMILeadingPhotonCandidateLen([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    return slc->reco.pfp[idx].shw.len;
  });

  const Var kNuMISubLeadingPhotonCandidateLen([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    return slc->reco.pfp[idx].shw.len;
  });

  const Var kPi0LeadingPhotonCandidateHitCompletenessBestmatch([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.bestmatch.hit_completeness) || std::isinf(slc->reco.pfp[idx].shw.truth.bestmatch.hit_completeness) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.bestmatch.hit_completeness;
  });

  const Var kPi0LeadingPhotonCandidateInFV([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    
    bool IsStartInFV = isInFV(slc->reco.pfp[idx].shw.start.x, slc->reco.pfp[idx].shw.start.y, slc->reco.pfp[idx].shw.start.z);
    bool IsEndInFV = isInFV(slc->reco.pfp[idx].shw.end.x, slc->reco.pfp[idx].shw.end.y, slc->reco.pfp[idx].shw.end.z);

    float IsShowerInFV = -5.f;
    if (IsStartInFV && IsEndInFV){
      IsShowerInFV = 1.f;
    }else if(IsStartInFV){
      IsShowerInFV = 0.f;
    }else{
      IsShowerInFV = -1.f;
    }
    return IsShowerInFV;
  });

  const Var kPi0LeadingPhotonCandidateEnergyCompletenessBestmatch([](const caf::SRSliceProxy* slc) -> float {
      int idx = kNuMILeadingPhotonCandidateIdx(slc);
      if(idx<0) return -5.f;
      if ( std::isnan(slc->reco.pfp[idx].shw.truth.bestmatch.energy_completeness) || std::isinf(slc->reco.pfp[idx].shw.truth.bestmatch.energy_completeness) ) return -5.f;
      return slc->reco.pfp[idx].shw.truth.bestmatch.energy_completeness;
  });


  const Var kPi0SubLeadingPhotonCandidateHitCompletenessBestmatch([](const caf::SRSliceProxy* slc) -> float {
      int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
      if(idx<0) return -5.f;
      if ( std::isnan(slc->reco.pfp[idx].shw.truth.bestmatch.hit_completeness) || std::isinf(slc->reco.pfp[idx].shw.truth.bestmatch.hit_completeness) ) return -5.f;
      return slc->reco.pfp[idx].shw.truth.bestmatch.hit_completeness;
  });

  const Var kPi0SubLeadingPhotonCandidateEnergyCompletenessBestmatch([](const caf::SRSliceProxy* slc) -> float {
      int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
      if(idx<0) return -5.f;
      if ( std::isnan(slc->reco.pfp[idx].shw.truth.bestmatch.energy_completeness) || std::isinf(slc->reco.pfp[idx].shw.truth.bestmatch.energy_completeness) ) return -5.f;
      return slc->reco.pfp[idx].shw.truth.bestmatch.energy_completeness;
  });

  const Var kPi0SubLeadingPhotonCandidateInFV([](const caf::SRSliceProxy* slc) -> float {
      int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
      if(idx<0) return -5.f;
      
      bool IsStartInFV = isInFV(slc->reco.pfp[idx].shw.start.x, slc->reco.pfp[idx].shw.start.y, slc->reco.pfp[idx].shw.start.z);
      bool IsEndInFV = isInFV(slc->reco.pfp[idx].shw.end.x, slc->reco.pfp[idx].shw.end.y, slc->reco.pfp[idx].shw.end.z);

      float IsShowerInFV = -5.f;
      if (IsStartInFV && IsEndInFV){
        IsShowerInFV = 1.f;
      }else if(IsStartInFV){
        IsShowerInFV = 0.f;
      }else{
        IsShowerInFV = -1.f;
      }
      return IsShowerInFV;
  });

  const Var kPi0LeadingPhotonCandidateBestmatchG4ID([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float bestmatchG4ID = slc->reco.pfp[idx].shw.truth.bestmatch.G4ID;
    if ( std::isnan(bestmatchG4ID) || std::isinf(bestmatchG4ID) ) return -5.f;
    return bestmatchG4ID;
  });

  const Var kPi0SubLeadingPhotonCandidateBestmatchG4ID([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float bestmatchG4ID = slc->reco.pfp[idx].shw.truth.bestmatch.G4ID;
    if ( std::isnan(bestmatchG4ID) || std::isinf(bestmatchG4ID) ) return -5.f;
    return bestmatchG4ID;
  });


  const Var kPi0LeadingPhotonCandidateG4ID([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float G4ID = slc->reco.pfp[idx].shw.truth.p.G4ID;
    if ( std::isnan(G4ID) || std::isinf(G4ID) ) return -5.f;
    return G4ID;
  });

  const Var kPi0SubLeadingPhotonCandidateG4ID([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float G4ID = slc->reco.pfp[idx].shw.truth.p.G4ID;
    if ( std::isnan(G4ID) || std::isinf(G4ID) ) return -5.f;
    return G4ID;
  });

  const Var kPi0LeadingPhotonCandidateBestplane_Energy([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float bestplane_energy = slc->reco.pfp[idx].shw.bestplane_energy;
    if ( std::isnan(bestplane_energy) || std::isinf(bestplane_energy) ) return -5.f;
    return bestplane_energy;
  });

  const Var kPi0SubLeadingPhotonCandidateBestplane_Energy([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float bestplane_energy = slc->reco.pfp[idx].shw.bestplane_energy;
    if ( std::isnan(bestplane_energy) || std::isinf(bestplane_energy) ) return -5.f;
    return bestplane_energy;
  });

  const Var kPi0LeadingPhotonCandidateBestplane_dEdx([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float bestplane_dEdx = slc->reco.pfp[idx].shw.bestplane_dEdx;
    if ( std::isnan(bestplane_dEdx) || std::isinf(bestplane_dEdx) ) return -5.f;
    return bestplane_dEdx;
  });

  const Var kPi0SubLeadingPhotonCandidateBestplane_dEdx([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float bestplane_dEdx = slc->reco.pfp[idx].shw.bestplane_dEdx;
    if ( std::isnan(bestplane_dEdx) || std::isinf(bestplane_dEdx) ) return -5.f;
    return bestplane_dEdx;
  });

  const Var kPi0LeadingPhotonCandidateIsContained([](const caf::SRSliceProxy* slc) -> float {
      int idx = kNuMILeadingPhotonCandidateIdx(slc);
      if(idx<0) return -5.f;

      bool IsStartCont = isContainedVol(slc->reco.pfp[idx].shw.start.x, slc->reco.pfp[idx].shw.start.y, slc->reco.pfp[idx].shw.start.z);
      bool IsEndCont = isContainedVol(slc->reco.pfp[idx].shw.end.x, slc->reco.pfp[idx].shw.end.y, slc->reco.pfp[idx].shw.end.z);

      float IsShowerCont = -5.f;
      if (IsStartCont && IsEndCont){
        IsShowerCont = 1.f;
      }else if(IsStartCont){
        IsShowerCont = 0.f;
      }else{
        IsShowerCont = -1.f;
      }
      return IsShowerCont;
  });

  const Var kPi0SubLeadingPhotonCandidateIsContained([](const caf::SRSliceProxy* slc) -> float {
      int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
      if(idx<0) return -5.f;      
      bool IsStartCont = isContainedVol(slc->reco.pfp[idx].shw.start.x, slc->reco.pfp[idx].shw.start.y, slc->reco.pfp[idx].shw.start.z);
      bool IsEndCont = isContainedVol(slc->reco.pfp[idx].shw.end.x, slc->reco.pfp[idx].shw.end.y, slc->reco.pfp[idx].shw.end.z);

      float IsShowerCont = -5.f;
      // 1 is fully contained, 0 is only start contained, -1 is not contained
      if (IsStartCont && IsEndCont){
        IsShowerCont = 1.f;
      }else if(IsStartCont){
        IsShowerCont = 0.f;
      }else{
        IsShowerCont = -1.f;
      }
      return IsShowerCont;
  });

  const Var kPi0LeadingPhotonCandidateTrackScore([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float trkscore = slc->reco.pfp[idx].trackScore;
    if ( std::isnan(trkscore) || std::isinf(trkscore) ) return -5.f;
    return trkscore;
  });

  const Var kPi0SubLeadingPhotonCandidateTrackScore([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float trkscore = slc->reco.pfp[idx].trackScore;
    if ( std::isnan(trkscore) || std::isinf(trkscore) ) return -5.f;
    return trkscore;
  });

  const Var kPi0LeadingPhotonCandidatePur([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float pur = slc->reco.pfp[idx].shw.truth.pur;
    if ( std::isnan(pur) || std::isinf(pur) ) return -5.f;
    return pur;
  });

  const Var kPi0SubLeadingPhotonCandidatePur([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float pur = slc->reco.pfp[idx].shw.truth.pur;
    if ( std::isnan(pur) || std::isinf(pur) ) return -5.f;
    return pur;
  });

  const Var kPi0LeadingPhotonCandidateEff([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float eff = slc->reco.pfp[idx].shw.truth.eff;
    if ( std::isnan(eff) || std::isinf(eff) ) return -5.f;
    return eff;
  });

  const Var kPi0SubLeadingPhotonCandidateEff([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float eff = slc->reco.pfp[idx].shw.truth.eff;
    if ( std::isnan(eff) || std::isinf(eff) ) return -5.f;
    return eff;
  });

const Var kPi0LeadingPhotonCandidateStartX([](const caf::SRSliceProxy* slc) -> float {
   int idx = kNuMILeadingPhotonCandidateIdx(slc);
   if(idx<0) return -99999.f;
   float start = slc->reco.pfp[idx].shw.start.x;
   if ( std::isnan(start) || std::isinf(start) ) return -99999.f;
   return start;
 });

const Var kPi0LeadingPhotonCandidateStartY([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMILeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  float start = slc->reco.pfp[idx].shw.start.y;
  if ( std::isnan(start) || std::isinf(start) ) return -99999.f;
  return start;
});

const Var kPi0LeadingPhotonCandidateStartZ([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMILeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  float start = slc->reco.pfp[idx].shw.start.z;
  if ( std::isnan(start) || std::isinf(start) ) return -99999.f;
  return start;
});

const Var kPi0SubLeadingPhotonCandidateStartX([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  float start = slc->reco.pfp[idx].shw.start.x;
  if ( std::isnan(start) || std::isinf(start) ) return -99999.f;
  return start;
});

const Var kPi0SubLeadingPhotonCandidateStartY([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
 if(idx<0) return -99999.f;
 float start = slc->reco.pfp[idx].shw.start.y;
 if ( std::isnan(start) || std::isinf(start) ) return -99999.f;
 return start;
});

const Var kPi0SubLeadingPhotonCandidateStartZ([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
 if(idx<0) return -99999.f;
 float start = slc->reco.pfp[idx].shw.start.z;
 if ( std::isnan(start) || std::isinf(start) ) return -99999.f;
 return start;
});

const Var kPi0LeadingPhotonCandidateTrueStartX([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMILeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.start.x) || std::isinf(slc->reco.pfp[idx].shw.truth.p.start.x) ) return -99999.f;
  return slc->reco.pfp[idx].shw.truth.p.start.x;
});

const Var kPi0LeadingPhotonCandidateTrueStartY([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMILeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.start.y) || std::isinf(slc->reco.pfp[idx].shw.truth.p.start.y) ) return -99999.f;
  return slc->reco.pfp[idx].shw.truth.p.start.y;
});

const Var kPi0LeadingPhotonCandidateTrueStartZ([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMILeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.start.z) || std::isinf(slc->reco.pfp[idx].shw.truth.p.start.z) ) return -99999.f;
  return slc->reco.pfp[idx].shw.truth.p.start.z;
});

const Var kPi0SubLeadingPhotonCandidateTrueStartX([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.start.x) || std::isinf(slc->reco.pfp[idx].shw.truth.p.start.x) ) return -99999.f;
  return slc->reco.pfp[idx].shw.truth.p.start.x;
});  

const Var kPi0SubLeadingPhotonCandidateTrueStartY([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.start.y) || std::isinf(slc->reco.pfp[idx].shw.truth.p.start.y) ) return -99999.f;
  return slc->reco.pfp[idx].shw.truth.p.start.y;
});  

const Var kPi0SubLeadingPhotonCandidateTrueStartZ([](const caf::SRSliceProxy* slc) -> float {
  int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
  if(idx<0) return -99999.f;
  if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.start.z) || std::isinf(slc->reco.pfp[idx].shw.truth.p.start.z) ) return -99999.f;
  return slc->reco.pfp[idx].shw.truth.p.start.z;
}); 

  const Var kPi0LeadingPhotonCandidateConversionGap([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -9999.f;
    float var = slc->reco.pfp[idx].shw.conversion_gap;
    if ( std::isnan(var) || std::isinf(var) ) return -9999.f;
    return var;
  });

  const Var kPi0SubLeadingPhotonCandidateConversionGap([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -9999.f;
    float var = slc->reco.pfp[idx].shw.conversion_gap;
    if ( std::isnan(var) || std::isinf(var) ) return -9999.f;
    return var;
  });

  const Var kBaryDeltaY([](const caf::SRSliceProxy *slc) -> double {
    return slc->barycenterFM.deltaY_Trigger;
  });

  const Var kBaryDeltaZ([](const caf::SRSliceProxy *slc) -> double {
    return slc->barycenterFM.deltaZ_Trigger;
  });

  const Var kBaryRadius([](const caf::SRSliceProxy *slc) -> double {
    return slc->barycenterFM.radius_Trigger;
  });

  const Var kBaryFlashFirstHit([](const caf::SRSliceProxy *slc) -> double {
    return slc->barycenterFM.flashFirstHit;
  });

  const Var kPi0LeadingPhotonCandidateCosmicDist([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float var = slc->reco.pfp[idx].shw.cosmicDist;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;
    return var;
  });

  const Var kPi0SubLeadingPhotonCandidateCosmicDist([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float var = slc->reco.pfp[idx].shw.cosmicDist;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;
    return var;
  });

  const Var kPi0LeadingPhotonCandidateShowerDensity([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float var = slc->reco.pfp[idx].shw.density;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;
    return var;
  });
  
  const Var kPi0SubLeadingPhotonCandidateShowerDensity([](const caf::SRSliceProxy* slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    float var = slc->reco.pfp[idx].shw.density;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;
    return var;
  });

  const Var kPi0LeadingPhotonCandidateShowerGenPX([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -99999.f;
    float var = slc->reco.pfp[idx].shw.truth.p.genp.x;
    if ( std::isnan(var) || std::isinf(var) ) return -99999.f;
    return var;
  });
  
  const Var kPi0SubLeadingPhotonCandidateShowerGenPX([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -99999.f;
    float var = slc->reco.pfp[idx].shw.truth.p.genp.x;
    if ( std::isnan(var) || std::isinf(var) ) return -99999.f;
    return var;
  });

  const Var kPi0LeadingPhotonCandidateShowerGenPY([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -99999.f;
    float var = slc->reco.pfp[idx].shw.truth.p.genp.y;
    if ( std::isnan(var) || std::isinf(var) ) return -99999.f;
    return var;
  });

  const Var kPi0SubLeadingPhotonCandidateShowerGenPY([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -99999.f;
    float var = slc->reco.pfp[idx].shw.truth.p.genp.y;
    if ( std::isnan(var) || std::isinf(var) ) return -99999.f;
    return var;
  });
  
  const Var kPi0LeadingPhotonCandidateShowerGenPZ([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -99999.f;
    float var = slc->reco.pfp[idx].shw.truth.p.genp.z;
    if ( std::isnan(var) || std::isinf(var) ) return -99999.f;
    return var;
  });

  const Var kPi0SubLeadingPhotonCandidateShowerGenPZ([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -99999.f;
    float var = slc->reco.pfp[idx].shw.truth.p.genp.z;
    if ( std::isnan(var) || std::isinf(var) ) return -99999.f;
    return var;
  });

  const Var kPi0LeadingPhotonCandidateTrueCryostat([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.cryostat) || std::isinf(slc->reco.pfp[idx].shw.truth.p.cryostat) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.cryostat;
  });

  const Var kPi0SubLeadingPhotonCandidateTrueCryostat([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.cryostat) || std::isinf(slc->reco.pfp[idx].shw.truth.p.cryostat) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.cryostat;
  });

  const Var kPi0LeadingPhotonCandidateTrueLength([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.length) || std::isinf(slc->reco.pfp[idx].shw.truth.p.length) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.length;
  });

  const Var kPi0SubLeadingPhotonCandidateTrueLength([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.length) || std::isinf(slc->reco.pfp[idx].shw.truth.p.length) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.length;
  });

  const Var kPi0LeadingPhotonCandidatePDG([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.pdg) || std::isinf(slc->reco.pfp[idx].shw.truth.p.pdg) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.pdg;
  });

  const Var kPi0SubLeadingPhotonCandidatePDG([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.truth.p.pdg) || std::isinf(slc->reco.pfp[idx].shw.truth.p.pdg) ) return -5.f;
    return slc->reco.pfp[idx].shw.truth.p.pdg;
  });

  const Var kPi0LeadingPhotonCandidateProtonChi2([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_proton) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_proton) ) return -5.f;
    return slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_proton;
  });

  const Var kPi0SubLeadingPhotonCandidateProtonChi2([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_proton) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_proton) ) return -5.f;
    return slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_proton;
  });

  const Var kPi0LeadingPhotonCandidateMuonChi2([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) ) return -5.f;
    return slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon;
  });

  const Var kPi0SubLeadingPhotonCandidateMuonChi2([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) ) return -5.f;
    return slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon;
  });

  const Var kPi0LeadingPhotonCandidatePionChi2([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) ) return -5.f;
    return slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion;
  });

  const Var kPi0SubLeadingPhotonCandidatePionChi2([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) ) return -5.f;
    return slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion;
  });

  const Var kPi0LeadingPhotonCandidateMuonPionChi2Diff([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) ) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) ) return -5.f;
    float diff = slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon - slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion;
    return diff;
  });

  const Var kPi0SubLeadingPhotonCandidateMuonPionChi2Diff([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon) ) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) || std::isinf(slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion) ) return -5.f;
    float diff = slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_muon - slc->reco.pfp[idx].trk.chi2pid[bestplane].chi2_pion;
    return diff;
  });

  const Var kPi0LeadingPhotonCandidateNHits([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].shw.plane[bestplane].nHits) || std::isinf(slc->reco.pfp[idx].shw.plane[bestplane].nHits) ) return -5.f;
    return slc->reco.pfp[idx].shw.plane[bestplane].nHits;
  });

  const Var kPi0SubLeadingPhotonCandidateNHits([](const caf::SRSliceProxy *slc) -> int {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].shw.plane[bestplane].nHits) || std::isinf(slc->reco.pfp[idx].shw.plane[bestplane].nHits) ) return -5.f;
    return slc->reco.pfp[idx].shw.plane[bestplane].nHits;
  });

  const Var kPi0LeadingPhotonCandidateSqrtEDen([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.bestplane_energy) || std::isinf(slc->reco.pfp[idx].shw.bestplane_energy) ) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.len) || std::isinf(slc->reco.pfp[idx].shw.len) ) return -5.f;
    float sqrtEDen = std::sqrt(slc->reco.pfp[idx].shw.bestplane_energy) / slc->reco.pfp[idx].shw.len;
    if ( std::isnan(sqrtEDen) || std::isinf(sqrtEDen) ) return -5.f;
    return sqrtEDen;
  });

  const Var kPi0SubLeadingPhotonCandidateSqrtEDen([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.bestplane_energy) || std::isinf(slc->reco.pfp[idx].shw.bestplane_energy) ) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.len) || std::isinf(slc->reco.pfp[idx].shw.len) ) return -5.f;
    float sqrtEDen = std::sqrt(slc->reco.pfp[idx].shw.bestplane_energy) / slc->reco.pfp[idx].shw.len;
    if ( std::isnan(sqrtEDen) || std::isinf(sqrtEDen) ) return -5.f;
    return sqrtEDen;
  });

  const Var kPi0LeadingPhotonCandidateHitDen([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMILeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].shw.plane[bestplane].nHits) || std::isinf(slc->reco.pfp[idx].shw.plane[bestplane].nHits) ) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.plane[bestplane].wirePitch) || std::isinf(slc->reco.pfp[idx].shw.plane[bestplane].wirePitch)) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.len) || std::isinf(slc->reco.pfp[idx].shw.len) ) return -5.f;
    float len = slc->reco.pfp[idx].shw.len, nhits = slc->reco.pfp[idx].shw.plane[bestplane].nHits, wirePitch = slc->reco.pfp[idx].shw.plane[bestplane].wirePitch;
    float wiresHit = len / wirePitch;
    if ( std::isnan(wiresHit) || std::isinf(wiresHit) ) return -5.f;
    float hitDen = nhits / wiresHit;
    if ( std::isnan(hitDen) || std::isinf(hitDen) ) return -5.f;
    return hitDen;
  });

  const Var kPi0SubLeadingPhotonCandidateHitDen([](const caf::SRSliceProxy *slc) -> float {
    int idx = kNuMISubLeadingPhotonCandidateIdx(slc);
    if(idx<0) return -5.f;
    int bestplane = slc->reco.pfp[idx].shw.bestplane;
    if ( std::isnan(slc->reco.pfp[idx].shw.plane[bestplane].nHits) || std::isinf(slc->reco.pfp[idx].shw.plane[bestplane].nHits) ) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.plane[bestplane].wirePitch) || std::isinf(slc->reco.pfp[idx].shw.plane[bestplane].wirePitch)) return -5.f;
    if ( std::isnan(slc->reco.pfp[idx].shw.len) || std::isinf(slc->reco.pfp[idx].shw.len) ) return -5.f;
    float len = slc->reco.pfp[idx].shw.len, nhits = slc->reco.pfp[idx].shw.plane[bestplane].nHits, wirePitch = slc->reco.pfp[idx].shw.plane[bestplane].wirePitch;
    float wiresHit = len / wirePitch;
    if ( std::isnan(wiresHit) || std::isinf(wiresHit) ) return -5.f;
    float hitDen = nhits / wiresHit;
    if ( std::isnan(hitDen) || std::isinf(hitDen) ) return -5.f;
    return hitDen;
  });

  const Var kIsClearCosmic([](const caf::SRSliceProxy *slc) -> float {
    return slc->is_clear_cosmic;
  });

}




