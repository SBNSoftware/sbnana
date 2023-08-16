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
        if ( (!Contained && trk.len > 100.) || (Contained && trk.len > 50. && Chi2Proton > 60. && Chi2Muon < 30.) ) {
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
    if ( slc->truth.index < 0 ) return p;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 && momentum > p ) {
        p = momentum;
      }
    }

    return p;
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
    if ( slc->truth.index < 0 ) return p;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 2212 && prim.contained && momentum > p ) {
        p = momentum;
      }
    }

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
    float p(-5.f);
    float cosThNuMI(-5.f);
    if ( slc->truth.index < 0 ) return cosThNuMI;

      double magNuMI = sqrt(315.120380*315.120380 + 33.644912*33.644912 + 733.632532*733.632532);
      TVector3 rFromNuMI(315.120380/magNuMI, 33.644912/magNuMI, 733.632532/magNuMI);

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 && momentum > p ) {
        p = momentum;

        TVector3 muDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
        muDir = muDir.Unit();
        cosThNuMI = TMath::Cos( muDir.Angle(rFromNuMI) );
      }
    }

    return cosThNuMI;
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
    float pM(-5.f);
    float pP(-5.f);
    float pX(0.f);
    float pY(0.f);
    float pZ(0.f);
    float mX(0.f);
    float mY(0.f);
    float mZ(0.f);
    if ( slc->truth.index < 0 ) return -5.f;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 && momentum > pM ) {
        pM = momentum;

        TVector3 muDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
        muDir = muDir.Unit();
        mX = muDir.X();
        mY = muDir.Y();
        mZ = muDir.Z();
      }
      if ( abs(prim.pdg) == 2212 && prim.contained && momentum > pP ) {
        pP = momentum;

        TVector3 pDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
        pDir = pDir.Unit();
        pX = pDir.X();
        pY = pDir.Y();
        pZ = pDir.Z();
      }
    }

    TVector3 muonDir(mX, mY, mZ);
    muonDir = muonDir.Unit();
    TVector3 protonDir(pX, pY, pZ);
    protonDir = protonDir.Unit();
    return TMath::Cos( muonDir.Angle(protonDir) );
  });

}
