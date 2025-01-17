#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnana/SBNAna/Vars/Pi0StudyVars.h"

#include "TVector3.h"

namespace ana {

  // Utility functions
  bool isInFV (double x, double y, double z)
  {
    if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

    return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
              ( x >  61.94 + 25 && x <  358.49 - 25 )) &&
            ( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
              ( z > -894.95 + 30 && z < 894.95 - 50 ) ));
  }

  bool isContainedVol (double x, double y, double z)
  {
    if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

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

  // Muon candidate
  const Var kNuMIMuonCandidateIdx([](const caf::SRSliceProxy* slc) -> int {
      float Longest(0);
      int PTrackInd(-1);

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
        const float Chi2Proton = trk.chi2pid[2].chi2_proton;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        //track len of muon
        if ( (!Contained && trk.len > 0.) || (Contained && trk.len > 0. && Chi2Proton > 60. && Chi2Muon < 30.) ) {
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
    int PTrackInd(-1);

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

  // MultiVar for the photon candidate indices
  const MultiVar kNuMIPhotonCandidateIdxs([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;

    int primaryInd = kNuMIMuonCandidateIdx(slc);
    int primaryProtonInd = kNuMIProtonCandidateIdx(slc);

    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); ++i_pfp){

      if ( i_pfp == (unsigned int)primaryInd || i_pfp == (unsigned int)primaryProtonInd ) {
        continue; // skip the particle which is the muon or leading proton candidate!
      }
      if ( !IsShowerlike(slc, i_pfp) ) { 
        continue; // skip things with track score > 0.45
      }
      auto const& shw = slc->reco.pfp.at(i_pfp).shw;

      // Check if shower fit even seems kind-of valid:
      if ( std::isnan(shw.start.x) || (shw.start.x > -5.5 && shw.start.x < -4.5) ||
           std::isnan(shw.len) || shw.len <= 0. ) continue;

      // if it meets this then we're not going to cut on it...
      if ( std::isnan(shw.plane[2].energy) || std::isinf(shw.plane[2].energy) || shw.plane[2].energy <= 0.0 ) continue;

      // and... if it meets then then we're not going to cut on it...
      if ( std::isnan(shw.conversion_gap) || std::isinf(shw.conversion_gap) || shw.conversion_gap <= 0. ) continue;

      // if we got here, then it should be the case that the fit seems valid and:
      // shwE > 0.040 GeV
      // trackScore < 0.45 (technically <= 0.45)
      // conversionGap > 5. cm

      rets.push_back( i_pfp );

    }

    // guess we're not cutting anything
    return rets;

  });

  const Var kNumberRecoShowers([](const caf::SRSliceProxy* slc) -> int {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc); //track score > 0.45
    if (std::isnan(photon_indices.size())) return 0;
    int numshowers = photon_indices.size();
    return numshowers;
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
      else p = trk.mcsP.fwdP_muon;
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

  const Var kNuMILeadingPhotonCandidateE([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()==0) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;

    return maxE;
  });

  const Var kNuMILeadingPhotonCandidateTrueE([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;

    if ( std::isnan(slc->reco.pfp[idxMaxE].shw.truth.p.genE) || std::isinf(slc->reco.pfp[idxMaxE].shw.truth.p.genE) ) return -5.f;
    float trueE = slc->reco.pfp[idxMaxE].shw.truth.p.genE;

    return trueE;
  });

  const Var kNuMISecondaryPhotonCandidateE([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;

    return scdy;
  });

  const Var kNuMISecondaryPhotonCandidateTrueE([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    
    if ( std::isnan(slc->reco.pfp[idxScdy].shw.truth.p.genE) || std::isinf(slc->reco.pfp[idxScdy].shw.truth.p.genE) ) return -5.f;
    float trueE = slc->reco.pfp[idxScdy].shw.truth.p.genE;
    return trueE;
  });

  const Var kNuMIPhotonCandidatesOpeningAngle([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( idxMaxE == idxScdy ) return -5.f;
    if ( maxE < 0. || scdy < 0. ) return -5.f;

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
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()==0) return -5.f;

    //Find length of most energetic shower
    float len = -5.;
    float maxE = -5.;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        len = slc->reco.pfp[idxI].shw.len;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;

    return len;
  });

  const Var kNuMISecondaryPhotonCandidateLen([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    float len;
    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        len = slc->reco.pfp[idxI].shw.len; 
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
        len = slc->reco.pfp[idxI].shw.len;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;

    return len;
  });

  const Var kPi0LeadingPhotonCandidateHitCompletenessBestmatch([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    //rec.slc.reco.pfp.shw.truth.bestmatch.hit_completeness
    if ( std::isnan(slc->reco.pfp[idxMaxE].shw.truth.bestmatch.hit_completeness) || std::isinf(slc->reco.pfp[idxMaxE].shw.truth.bestmatch.hit_completeness) ) return -5.f;
    float maxPhotonCompleteness = slc->reco.pfp[idxMaxE].shw.truth.bestmatch.hit_completeness;
    return maxPhotonCompleteness;
  });

  const Var kPi0LeadingPhotonCandidateInFV([](const caf::SRSliceProxy* slc) -> float {
      std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
      if(photon_indices.size()<=1) return -5.f;

      // Find 2 most energetic:
      unsigned int idxMaxE = 0;
      float maxE = -5.;
      unsigned int idxScdy = 0;
      float scdy = -5.;
      for ( auto const& photon_idx : photon_indices ) {
        unsigned int idxI = (unsigned int)std::lround(photon_idx);
        if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
          idxScdy = idxMaxE;
          scdy = maxE;
          idxMaxE = idxI;
          maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
        else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
          idxScdy = idxI;
          scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
      }

      if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5;
      if (std::isnan(slc->reco.pfp[idxMaxE].shw.start.x) || std::isnan(slc->reco.pfp[idxMaxE].shw.start.y) || std::isnan(slc->reco.pfp[idxMaxE].shw.start.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxMaxE].shw.start.x) || std::isinf(slc->reco.pfp[idxMaxE].shw.start.y) || std::isinf(slc->reco.pfp[idxMaxE].shw.start.z)) return -5.f;

      if (std::isnan(slc->reco.pfp[idxMaxE].shw.end.x) || std::isnan(slc->reco.pfp[idxMaxE].shw.end.y) || std::isnan(slc->reco.pfp[idxMaxE].shw.end.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxMaxE].shw.end.x) || std::isinf(slc->reco.pfp[idxMaxE].shw.end.y) || std::isinf(slc->reco.pfp[idxMaxE].shw.end.z)) return -5.f;
      
      bool IsStartInFV = isInFV(slc->reco.pfp[idxMaxE].shw.start.x, slc->reco.pfp[idxMaxE].shw.start.y, slc->reco.pfp[idxMaxE].shw.start.z);
      bool IsEndInFV = isInFV(slc->reco.pfp[idxMaxE].shw.end.x, slc->reco.pfp[idxMaxE].shw.end.y, slc->reco.pfp[idxMaxE].shw.end.z);

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
      std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
      if(photon_indices.size()<=1) return -5.f;

      // Find 2 most energetic:
      unsigned int idxMaxE = 0;
      float maxE = -5.;
      unsigned int idxScdy = 0;
      float scdy = -5.;
      for ( auto const& photon_idx : photon_indices ) {
        unsigned int idxI = (unsigned int)std::lround(photon_idx);
        if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
          idxScdy = idxMaxE;
          scdy = maxE;
          idxMaxE = idxI;
          maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
        else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
          idxScdy = idxI;
          scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
      }

      if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
       //rec.slc.reco.pfp.shw.truth.bestmatch.hit_completeness
      if ( std::isnan(slc->reco.pfp[idxMaxE].shw.truth.bestmatch.energy_completeness) || std::isinf(slc->reco.pfp[idxMaxE].shw.truth.bestmatch.energy_completeness) ) return -5.f;
      float maxPhotonCompleteness = slc->reco.pfp[idxMaxE].shw.truth.bestmatch.energy_completeness;

      return maxPhotonCompleteness;
  });


  const Var kPi0SecondaryPhotonCandidateHitCompletenessBestmatch([](const caf::SRSliceProxy* slc) -> float {
      std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
      if(photon_indices.size()<=1) return -5.f;

      // Find 2 most energetic:
      unsigned int idxMaxE = 0;
      float maxE = -5.;
      unsigned int idxScdy = 0;
      float scdy = -5.;
      for ( auto const& photon_idx : photon_indices ) {
        unsigned int idxI = (unsigned int)std::lround(photon_idx);
        if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
          idxScdy = idxMaxE;
          scdy = maxE;
          idxMaxE = idxI;
          maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
        else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
          idxScdy = idxI;
          scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
      }

      if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;

      //rec.slc.reco.pfp.shw.truth.bestmatch.hit_completeness
      if ( std::isnan(slc->reco.pfp[idxScdy].shw.truth.bestmatch.hit_completeness) || std::isinf(slc->reco.pfp[idxScdy].shw.truth.bestmatch.hit_completeness) ) return -5.f;
      float maxPhotonCompleteness = slc->reco.pfp[idxScdy].shw.truth.bestmatch.hit_completeness;

      return maxPhotonCompleteness;
  });

  const Var kPi0SecondaryPhotonCandidateEnergyCompletenessBestmatch([](const caf::SRSliceProxy* slc) -> float {
        std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
        if(photon_indices.size()<=1) return -5.f;

        // Find 2 most energetic:
        unsigned int idxMaxE = 0;
        float maxE = -5.;
        unsigned int idxScdy = 0;
        float scdy = -5.;
        for ( auto const& photon_idx : photon_indices ) {
          unsigned int idxI = (unsigned int)std::lround(photon_idx);
          if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
            idxScdy = idxMaxE;
            scdy = maxE;
            idxMaxE = idxI;
            maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
          }
          else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
            idxScdy = idxI;
            scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
          }
        }

        if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;

         //rec.slc.reco.pfp.shw.truth.bestmatch.hit_completeness
        if ( std::isnan(slc->reco.pfp[idxScdy].shw.truth.bestmatch.energy_completeness) || std::isinf(slc->reco.pfp[idxScdy].shw.truth.bestmatch.energy_completeness) ) return -5.f;
        float maxPhotonCompleteness = slc->reco.pfp[idxScdy].shw.truth.bestmatch.energy_completeness;

        return maxPhotonCompleteness;
  });

  const Var kPi0SecondaryPhotonCandidateInFV([](const caf::SRSliceProxy* slc) -> float {
      std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
      if(photon_indices.size()<=1) return -5.f;

      // Find 2 most energetic:
      unsigned int idxMaxE = 0;
      float maxE = -5.;
      unsigned int idxScdy = 0;
      float scdy = -5.;
      for ( auto const& photon_idx : photon_indices ) {
        unsigned int idxI = (unsigned int)std::lround(photon_idx);
        if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
          idxScdy = idxMaxE;
          scdy = maxE;
          idxMaxE = idxI;
          maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
        else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
          idxScdy = idxI;
          scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
      }

      if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5;
      if (std::isnan(slc->reco.pfp[idxScdy].shw.start.x) || std::isnan(slc->reco.pfp[idxScdy].shw.start.y) || std::isnan(slc->reco.pfp[idxScdy].shw.start.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxScdy].shw.start.x) || std::isinf(slc->reco.pfp[idxScdy].shw.start.y) || std::isinf(slc->reco.pfp[idxScdy].shw.start.z)) return -5.f;

      if (std::isnan(slc->reco.pfp[idxScdy].shw.end.x) || std::isnan(slc->reco.pfp[idxScdy].shw.end.y) || std::isnan(slc->reco.pfp[idxScdy].shw.end.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxScdy].shw.end.x) || std::isinf(slc->reco.pfp[idxScdy].shw.end.y) || std::isinf(slc->reco.pfp[idxScdy].shw.end.z)) return -5.f;
      
      bool IsStartInFV = isInFV(slc->reco.pfp[idxScdy].shw.start.x, slc->reco.pfp[idxScdy].shw.start.y, slc->reco.pfp[idxScdy].shw.start.z);
      bool IsEndInFV = isInFV(slc->reco.pfp[idxScdy].shw.end.x, slc->reco.pfp[idxScdy].shw.end.y, slc->reco.pfp[idxScdy].shw.end.z);

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
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float bestmatchG4ID = slc->reco.pfp[idxMaxE].shw.truth.bestmatch.G4ID;
    if ( std::isnan(bestmatchG4ID) || std::isinf(bestmatchG4ID) ) return -5.f;

    return bestmatchG4ID;
  });

  const Var kPi0SecondaryPhotonCandidateBestmatchG4ID([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float bestmatchG4ID = slc->reco.pfp[idxScdy].shw.truth.bestmatch.G4ID;
    if ( std::isnan(bestmatchG4ID) || std::isinf(bestmatchG4ID) ) return -5.f;

    return bestmatchG4ID;
  });


  const Var kPi0LeadingPhotonCandidateG4ID([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float G4ID = slc->reco.pfp[idxMaxE].shw.truth.p.G4ID;
    if ( std::isnan(G4ID) || std::isinf(G4ID) ) return -5.f;

    return G4ID;
  });

  const Var kPi0SecondaryPhotonCandidateG4ID([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float G4ID = slc->reco.pfp[idxScdy].shw.truth.p.G4ID;
    if ( std::isnan(G4ID) || std::isinf(G4ID) ) return -5.f;

    return G4ID;
  });

  const Var kPi0LeadingPhotonCandidateBestplane_Energy([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float bestplane_energy = slc->reco.pfp[idxMaxE].shw.bestplane_energy;
    if ( std::isnan(bestplane_energy) || std::isinf(bestplane_energy) ) return -5.f;

    return bestplane_energy;
  });

  const Var kPi0SecondaryPhotonCandidateBestplane_Energy([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float bestplane_energy = slc->reco.pfp[idxScdy].shw.bestplane_energy;
    if ( std::isnan(bestplane_energy) || std::isinf(bestplane_energy) ) return -5.f;

    return bestplane_energy;
  });

  const Var kPi0LeadingPhotonCandidateBestplane_dEdx([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float bestplane_dEdx = slc->reco.pfp[idxMaxE].shw.bestplane_dEdx;
    if ( std::isnan(bestplane_dEdx) || std::isinf(bestplane_dEdx) ) return -5.f;

    return bestplane_dEdx;
  });

  const Var kPi0SecondaryPhotonCandidateBestplane_dEdx([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float bestplane_dEdx = slc->reco.pfp[idxScdy].shw.bestplane_dEdx;
    if ( std::isnan(bestplane_dEdx) || std::isinf(bestplane_dEdx) ) return -5.f;

    return bestplane_dEdx;
  });

  const Var kPi0LeadingPhotonCandidateIsContained([](const caf::SRSliceProxy* slc) -> float {
      std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
      if(photon_indices.size()<=1) return -5.f;

      // Find 2 most energetic:
      unsigned int idxMaxE = 0;
      float maxE = -5.;
      unsigned int idxScdy = 0;
      float scdy = -5.;
      for ( auto const& photon_idx : photon_indices ) {
        unsigned int idxI = (unsigned int)std::lround(photon_idx);
        if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
          idxScdy = idxMaxE;
          scdy = maxE;
          idxMaxE = idxI;
          maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
        else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
          idxScdy = idxI;
          scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
      }

      if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5;
      if (std::isnan(slc->reco.pfp[idxMaxE].shw.start.x) || std::isnan(slc->reco.pfp[idxMaxE].shw.start.y) || std::isnan(slc->reco.pfp[idxMaxE].shw.start.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxMaxE].shw.start.x) || std::isinf(slc->reco.pfp[idxMaxE].shw.start.y) || std::isinf(slc->reco.pfp[idxMaxE].shw.start.z)) return -5.f;

      if (std::isnan(slc->reco.pfp[idxMaxE].shw.end.x) || std::isnan(slc->reco.pfp[idxMaxE].shw.end.y) || std::isnan(slc->reco.pfp[idxMaxE].shw.end.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxMaxE].shw.end.x) || std::isinf(slc->reco.pfp[idxMaxE].shw.end.y) || std::isinf(slc->reco.pfp[idxMaxE].shw.end.z)) return -5.f;
      
      bool IsStartCont = isContainedVol(slc->reco.pfp[idxMaxE].shw.start.x, slc->reco.pfp[idxMaxE].shw.start.y, slc->reco.pfp[idxMaxE].shw.start.z);
      bool IsEndCont = isContainedVol(slc->reco.pfp[idxMaxE].shw.end.x, slc->reco.pfp[idxMaxE].shw.end.y, slc->reco.pfp[idxMaxE].shw.end.z);

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

  const Var kPi0SecondaryPhotonCandidateIsContained([](const caf::SRSliceProxy* slc) -> float {
      std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
      if(photon_indices.size()<=1) return -5.f;

      // Find 2 most energetic:
      unsigned int idxMaxE = 0;
      float maxE = -5.;
      unsigned int idxScdy = 0;
      float scdy = -5.;
      for ( auto const& photon_idx : photon_indices ) {
        unsigned int idxI = (unsigned int)std::lround(photon_idx);
        if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
          idxScdy = idxMaxE;
          scdy = maxE;
          idxMaxE = idxI;
          maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
        else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
          idxScdy = idxI;
          scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
        }
      }

      if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5;
      if (std::isnan(slc->reco.pfp[idxScdy].shw.start.x) || std::isnan(slc->reco.pfp[idxScdy].shw.start.y) || std::isnan(slc->reco.pfp[idxScdy].shw.start.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxScdy].shw.start.x) || std::isinf(slc->reco.pfp[idxScdy].shw.start.y) || std::isinf(slc->reco.pfp[idxScdy].shw.start.z)) return -5.f;

      if (std::isnan(slc->reco.pfp[idxScdy].shw.end.x) || std::isnan(slc->reco.pfp[idxScdy].shw.end.y) || std::isnan(slc->reco.pfp[idxScdy].shw.end.z)) return -5.f;
      if (std::isinf(slc->reco.pfp[idxScdy].shw.end.x) || std::isinf(slc->reco.pfp[idxScdy].shw.end.y) || std::isinf(slc->reco.pfp[idxScdy].shw.end.z)) return -5.f;
      
      bool IsStartCont = isContainedVol(slc->reco.pfp[idxScdy].shw.start.x, slc->reco.pfp[idxScdy].shw.start.y, slc->reco.pfp[idxScdy].shw.start.z);
      bool IsEndCont = isContainedVol(slc->reco.pfp[idxScdy].shw.end.x, slc->reco.pfp[idxScdy].shw.end.y, slc->reco.pfp[idxScdy].shw.end.z);

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


  const Var kPi0LeadingPhotonCandidateTrackScore([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float trkscore = slc->reco.pfp[idxMaxE].trackScore;
    if ( std::isnan(trkscore) || std::isinf(trkscore) ) return -5.f;

    return trkscore;
  });

  const Var kPi0SecondaryPhotonCandidateTrackScore([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float trkscore = slc->reco.pfp[idxScdy].trackScore;
    if ( std::isnan(trkscore) || std::isinf(trkscore) ) return -5.f;

    return trkscore;
  });

  const Var kPi0LeadingPhotonCandidatePur([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float pur = slc->reco.pfp[idxMaxE].shw.truth.pur;
    if ( std::isnan(pur) || std::isinf(pur) ) return -5.f;

    return pur;
  });

  const Var kPi0SecondaryPhotonCandidatePur([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float pur = slc->reco.pfp[idxScdy].shw.truth.pur;
    if ( std::isnan(pur) || std::isinf(pur) ) return -5.f;

    return pur;
  });


  const Var kPi0LeadingPhotonCandidateEff([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float eff = slc->reco.pfp[idxMaxE].shw.truth.eff;
    if ( std::isnan(eff) || std::isinf(eff) ) return -5.f;

    return eff;
  });

  const Var kPi0SecondaryPhotonCandidateEff([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float eff = slc->reco.pfp[idxScdy].shw.truth.eff;
    if ( std::isnan(eff) || std::isinf(eff) ) return -5.f;

    return eff;
  });


const Var kPi0LeadingPhotonCandidateStartX([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float start = slc->reco.pfp[idxMaxE].shw.start.x;
    if ( std::isnan(start) || std::isinf(start) ) return -99999.f;

    return start;
  });

const Var kPi0LeadingPhotonCandidateStartY([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    float start = slc->reco.pfp[idxMaxE].shw.start.y;
    if ( std::isnan(start) || std::isinf(start) ) return -5.f;

    return start;
  });

const Var kPi0LeadingPhotonCandidateStartZ([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    float start = slc->reco.pfp[idxMaxE].shw.start.z;
    if ( std::isnan(start) || std::isinf(start) ) return -5.f;

    return start;
  });

const Var kPi0SecondaryPhotonCandidateStartX([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    float start = slc->reco.pfp[idxScdy].shw.start.x;
    if ( std::isnan(start) || std::isinf(start) ) return -5.f;

    return start;
  });

const Var kPi0SecondaryPhotonCandidateStartY([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    float start = slc->reco.pfp[idxScdy].shw.start.y;
    if ( std::isnan(start) || std::isinf(start) ) return -5.f;

    return start;
  });

const Var kPi0SecondaryPhotonCandidateStartZ([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    float start = slc->reco.pfp[idxScdy].shw.start.z;
    if ( std::isnan(start) || std::isinf(start) ) return -5.f;

    return start;
  });


const Var kPi0LeadingPhotonCandidateTrueStartX([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    if ( std::isnan(slc->reco.pfp[idxMaxE].shw.truth.p.start.x) || std::isinf(slc->reco.pfp[idxMaxE].shw.truth.p.start.x) ) return -99999.f;
    float start = slc->reco.pfp[idxMaxE].shw.truth.p.start.x;

    return start;
  });

const Var kPi0LeadingPhotonCandidateTrueStartY([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    if ( std::isnan(slc->reco.pfp[idxMaxE].shw.truth.p.start.y) || std::isinf(slc->reco.pfp[idxMaxE].shw.truth.p.start.y) ) return -99999.f;
    float start = slc->reco.pfp[idxMaxE].shw.truth.p.start.y;

    return start;
  });

const Var kPi0LeadingPhotonCandidateTrueStartZ([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    if ( std::isnan(slc->reco.pfp[idxMaxE].shw.truth.p.start.z) || std::isinf(slc->reco.pfp[idxMaxE].shw.truth.p.start.z) ) return -99999.f;
    float start = slc->reco.pfp[idxMaxE].shw.truth.p.start.z;

    return start;
  });

const Var kPi0SecondaryPhotonCandidateTrueStartX([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    if ( std::isnan(slc->reco.pfp[idxScdy].shw.truth.p.start.x) || std::isinf(slc->reco.pfp[idxScdy].shw.truth.p.start.x) ) return -99999.f;
    float start = slc->reco.pfp[idxScdy].shw.truth.p.start.x;

    return start;
  });

const Var kPi0SecondaryPhotonCandidateTrueStartY([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    if ( std::isnan(slc->reco.pfp[idxScdy].shw.truth.p.start.y) || std::isinf(slc->reco.pfp[idxScdy].shw.truth.p.start.y) ) return -99999.f;
    float start = slc->reco.pfp[idxScdy].shw.truth.p.start.y;

    return start;
  });

const Var kPi0SecondaryPhotonCandidateTrueStartZ([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -99999.f;
    if ( std::isnan(slc->reco.pfp[idxScdy].shw.truth.p.start.z) || std::isinf(slc->reco.pfp[idxScdy].shw.truth.p.start.z) ) return -99999.f;
    float start = slc->reco.pfp[idxScdy].shw.truth.p.start.z;

    return start;
  });

  const Var kPi0LeadingPhotonCandidateConversionGap([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -9999.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float var = slc->reco.pfp[idxMaxE].shw.conversion_gap;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;

    return var;
  });

  const Var kPi0SecondaryPhotonCandidateConversionGap([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float var = slc->reco.pfp[idxScdy].shw.conversion_gap;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;

    return var;
  });

  const Var kBaryDeltaY([](const caf::SRSliceProxy *slc) -> double{
    return slc->barycenterFM.deltaY_Trigger;
  });

  const Var kBaryDeltaZ([](const caf::SRSliceProxy *slc) -> double{
    return slc->barycenterFM.deltaZ_Trigger;
  });

  const Var kBaryRadius([](const caf::SRSliceProxy *slc) -> double{
    return slc->barycenterFM.radius_Trigger;
  });

  const Var kBaryFlashFirstHit([](const caf::SRSliceProxy *slc) -> double{
    return slc->barycenterFM.flashFirstHit;
  });

  const Var kPi0LeadingPhotonCandidateCosmicDist([](const caf::SRSliceProxy *slc) -> float{
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float var = slc->reco.pfp[idxMaxE].shw.cosmicDist;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;

    return var;
  });

  const Var kPi0SecondaryPhotonCandidateCosmicDist([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float var = slc->reco.pfp[idxScdy].shw.cosmicDist;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;

    return var;
  });

   const Var kPi0LeadingPhotonCandidateShowerDensity([](const caf::SRSliceProxy *slc) -> float{
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }
    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float var = slc->reco.pfp[idxMaxE].shw.density;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;

    return var;
  });

  const Var kPi0SecondaryPhotonCandidateShowerDensity([](const caf::SRSliceProxy* slc) -> float {
    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()<=1) return -5.f;

    // Find 2 most energetic:
    unsigned int idxMaxE = 0;
    float maxE = -5.;
    unsigned int idxScdy = 0;
    float scdy = -5.;
    for ( auto const& photon_idx : photon_indices ) {
      unsigned int idxI = (unsigned int)std::lround(photon_idx);
      if ( slc->reco.pfp[idxI].shw.plane[2].energy > maxE ) {
        idxScdy = idxMaxE;
        scdy = maxE;
        idxMaxE = idxI;
        maxE = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
      else if ( slc->reco.pfp[idxI].shw.plane[2].energy > scdy ) {
        idxScdy = idxI;
        scdy = slc->reco.pfp[idxI].shw.plane[2].energy;
      }
    }

    if ( photon_indices.size()>=2 && idxMaxE == idxScdy ) return -5.f;
    float var = slc->reco.pfp[idxScdy].shw.density;
    if ( std::isnan(var) || std::isinf(var) ) return -5.f;

    return var;
  });


  
}




