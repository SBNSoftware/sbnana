#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"

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

    return (( ( x < -61.94 - 10. && x > -358.49 + 10. ) ||
              ( x >  61.94 + 10. && x <  358.49 - 10. )) &&
            ( ( y > -181.86 + 10. && y < 134.96 - 10. ) &&
              ( z > -894.95 + 10. && z < 894.95 - 10. ) ));
  }

  bool IsValidTrkIdx( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
    return slice->reco.npfp > idxTrk;
  }

  bool IsTracklikeTrack( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
    return (!std::isnan(slice->reco.pfp.at(idxTrk).trackScore) && slice->reco.pfp.at(idxTrk).trackScore > 0.45);
  }

  bool IsShowerlike( const caf::SRSliceProxy* slice, const unsigned int idxShw ) {
    return (!std::isnan(slice->reco.pfp.at(idxShw).trackScore) && slice->reco.pfp.at(idxShw).trackScore > 0. && slice->reco.pfp.at(idxShw).trackScore <= 0.45 );
  }

  bool IsPrimaryPFP( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
    return slice->reco.pfp.at(idxTrk).parent_is_primary;
  }

  // Emulated trigger time
  const SpillVar kNuMISpillTriggerTime ( [](const caf::SRSpillProxy *sr) -> double {
    double triggerTime = 0.;

    double foundTriggerTime = sr->hdr.triggerinfo.trigger_within_gate;
    if ( !std::isnan(foundTriggerTime) && !std::isinf(foundTriggerTime) && foundTriggerTime < -15. ) triggerTime = -15.;
    else if ( std::isnan(foundTriggerTime) || std::isinf(foundTriggerTime) ) triggerTime = -16.;
    else if ( foundTriggerTime > 30. ) triggerTime = 30.;
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
        if ( (!Contained && trk.len > 50.) || (Contained && trk.len > 50. && Chi2Proton > 60. && Chi2Muon < 30.) ) {
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

  // MultiVar for the proton candidate indices
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
      if ( std::isnan(shw.plane[2].energy) || std::isinf(shw.plane[2].energy) || shw.plane[2].energy <= 0.04 ) continue;

      // and... if it meets then then we're not going to cut on it...
      if ( std::isnan(shw.conversion_gap) || std::isinf(shw.conversion_gap) || shw.conversion_gap <= 5. ) continue;

      // if we got here, then it should be the case that the fit seems valid and:
      // shwE > 0.040 GeV
      // trackScore < 0.45 (technically <= 0.45)
      // conversionGap > 5. cm

      rets.push_back( i_pfp );

    }

    // guess we're not cutting anything
    return rets;

  });

  // Reco muon momentum
  const Var kNuMIMuonCandidateRecoP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kNuMIMuonCandidateIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
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


}
