#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  const Var MuonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
    if(primTrackIndices.size()==0){
      return -999.;
    }
    else{
      float Longest(0);
      int PTrackInd(-1);
      for(const auto& trkIdx: primTrackIndices){
        const auto& pfp = slc->reco.pfp.at(trkIdx);
        const auto& trk = pfp.trk;

        if(trk.bestplane == -1) continue;
        if(isnan(trk.start.x)) continue;
        // First we calculate the distance of each track to the slice vertex.
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

        const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
        const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

        const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
        {
          Longest = trk.len;
          PTrackInd = trkIdx;
        }

      }

      return PTrackInd;

    }
  });

  const Var MuonTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      return trk.len;
    }
    else{
      return -999.;
    }

  });
  const Var MuonTrackLengthMatchMuon([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(abs(this_pdg)!=13) return -999.;
      return MuonTrackLength(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackLengthMatchPionPlus([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=+211) return -999.;
      return MuonTrackLength(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackLengthMatchPionMinus([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=-211) return -999.;
      return MuonTrackLength(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackLengthMatchProton([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=2212) return -999.;
      return MuonTrackLength(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackLengthMatchElse([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(abs(this_pdg)==13 || abs(this_pdg)==211 || this_pdg==2212) return -999;
      return MuonTrackLength(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        if(trk.len>50.) return trk.rangeP.p_muon;
        else return -999.;
      }
      else{
        if(trk.len>100.) return trk.mcsP.fwdP_muon;
        else return -999.;
      }
    }
    else{
      return -999.;
    }

  });
  const Var MuonTrackPMatchMuon([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(abs(this_pdg)!=13) return -999.;
      return MuonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackPMatchPionPlus([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=+211) return -999.;
      return MuonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackPMatchPionMinus([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=-211) return -999.;
      return MuonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackPMatchProton([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=2212) return -999.;
      return MuonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackPMatchElse([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(abs(this_pdg)==13 || abs(this_pdg)==211 || this_pdg==2212) return -999;
      return MuonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& vtx = slc->vertex;

      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      //static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      //static const TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);

      //TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;
      TVector3 v3_numi_to_vtx = NuDirection_NuMI;

      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);

    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackNuMIToVtxCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& vtx = slc->vertex;

      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);

      TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;

      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);

    }
    else{
      return -999.;
    }
  });

  const Var ProtonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
    if(primTrackIndices.size()==0){
      return -999.;
    }
    else{
      float Longest(0);
      int PTrackInd(-1);
      int muonTrackIndex = MuonTrackIndex(slc);
      for(const auto& trkIdx: primTrackIndices){
        const auto& pfp = slc->reco.pfp.at(trkIdx);
        const auto& trk = pfp.trk;

        if(trkIdx==muonTrackIndex) continue;
        if(trk.bestplane == -1) continue;

        // First we calculate the distance of each track to the slice vertex.
        if(isnan(trk.start.x)) continue;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
/*
        const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
        const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
*/
        // pid from collection only
        const float Chi2Proton = trk.chi2pid[2].chi2_proton;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

        float angle = -5.0;
        if ( muonTrackIndex >= 0 ) {
          const unsigned int idxPrim = (unsigned int)muonTrackIndex;
          TVector3 muDir( slc->reco.pfp[idxPrim].trk.dir.x, slc->reco.pfp[idxPrim].trk.dir.y, slc->reco.pfp[idxPrim].trk.dir.z );
          TVector3 pDir( slc->reco.pfp[trkIdx].trk.dir.x, slc->reco.pfp[trkIdx].trk.dir.y, slc->reco.pfp[trkIdx].trk.dir.z );
          angle = TMath::Cos(muDir.Angle(pDir));
        }

        if ( AtSlice && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 && trk.len > Longest ) {
          Longest = trk.len;
          PTrackInd = trkIdx;
        }

      }

      return PTrackInd;

    }
  });


  const Var ProtonTrackP([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained) return trk.rangeP.p_proton;
      else return -999.;
    }
    else{
      return -999.;
    }

  });
  const Var ProtonTrackPMatchMuon([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(abs(this_pdg)!=13) return -999.;
      return ProtonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackPMatchPionPlus([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=+211) return -999.;
      return ProtonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackPMatchPionMinus([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=-211) return -999.;
      return ProtonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackPMatchProton([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(this_pdg!=2212) return -999.;
      return ProtonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackPMatchElse([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int this_pdg = trk.truth.p.pdg;
      if(abs(this_pdg)==13 || abs(this_pdg)==211 || this_pdg==2212) return -999.;
      return ProtonTrackP(slc);
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& vtx = slc->vertex;

      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);

      TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;

      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);

    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackNuMIToVtxCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& vtx = slc->vertex;

      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      //static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      //static const TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);

      //TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;
      TVector3 v3_numi_to_vtx = NuDirection_NuMI;

      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);

    }
    else{
      return -999.;
    }
  });

  namespace Aux{

    const Var RelaxedMuonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
      vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
      if(primTrackIndices.size()==0){
        return -999.;
      }
      else{

        // Requiring good proton track 
        int protonIdx = -1;
        double protonMaxLength = -1;
        for(const auto& trkIdx: primTrackIndices){
          const auto& pfp = slc->reco.pfp.at(trkIdx);
          const auto& trk = pfp.trk;

          if(trk.bestplane == -1) continue;

          // First we calculate the distance of each track to the slice vertex.
          if(isnan(trk.start.x)) continue;
          const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                         slc->vertex.y - trk.start.y,
                                         slc->vertex.z - trk.start.z);

          // We require that the distance of the track from the slice is less than
          // 10 cm and that the parent of the track has been marked as the primary.
          const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

          const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
          const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

          const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

          if ( AtSlice && Contained && Chi2Proton <= 100 && Chi2Muon >= 30 && trk.len > protonMaxLength ) {
            protonIdx = trkIdx;
            protonMaxLength = trk.len;
          }

        }

        int muonIdx = -1;
        if(protonIdx>=0){

          // then find the longest track

          double muonMaxLength = -1;
          for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
            const auto& pfp = slc->reco.pfp.at(i_pfp);
            const auto& trk = pfp.trk;

            if(i_pfp==(unsigned int)protonIdx) continue;
            if(trk.bestplane == -1) continue;

            // First we calculate the distance of each track to the slice vertex.
            if(isnan(trk.start.x)) continue;
            const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                           slc->vertex.y - trk.start.y,
                                           slc->vertex.z - trk.start.z);

            // We require that the distance of the track from the slice is less than
            // 10 cm and that the parent of the track has been marked as the primary.
            const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

            if ( AtSlice && trk.len > muonMaxLength ) {
              muonIdx = i_pfp;
              muonMaxLength = trk.len;
            }

          }

        }

        return muonIdx;

      }
    });
    const Var RelaxedMuonTrackTrackScore([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trackScore;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained) return trk.chi2pid[trk.bestplane].chi2_muon;
        else return -999.;

      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained) return trk.chi2pid[trk.bestplane].chi2_proton;
        else return -999.;

      }
      else{
        return -999.;
      }
    });

    const Var RelaxedProtonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
      vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
      if(primTrackIndices.size()==0){
        return -999.;
      }
      else{

        // Requiring good muon
        int muonTrackIndex = MuonTrackIndex(slc);
        int protonIdx = -1;
        bool HasPionCand = false;
        if(muonTrackIndex>=0){

          double protonMaxLength = -1.;

          // find the longest track
          // but the slice should not have another pion candidate

          for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
            const auto& pfp = slc->reco.pfp.at(i_pfp);
            const auto& trk = pfp.trk;

            if(i_pfp==(unsigned int )muonTrackIndex) continue;
            if(trk.bestplane == -1) continue;
            if(isnan(trk.start.x)) continue;
            // First we calculate the distance of each track to the slice vertex.
            const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                           slc->vertex.y - trk.start.y,
                                           slc->vertex.z - trk.start.z);

            // We require that the distance of the track from the slice is less than
            // 10 cm and that the parent of the track has been marked as the primary.
            const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

            const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
            const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;
            const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
            const bool MaybeTrackExiting = ( !Contained && trk.len > 100);
            const bool MaybeTrackContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
            // check pion cand other than the good muon we found
            if ( AtSlice && ( MaybeTrackExiting || MaybeTrackContained ) ){
              HasPionCand = true;
            }
            else{

              if( AtSlice && trk.len > protonMaxLength ){
                protonMaxLength = trk.len;
                protonIdx = i_pfp;
              }

            }

          }

        }

        if(HasPionCand) return -1;
        return protonIdx;

      }
    });
    const Var RelaxedProtonTrackTrackScore([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trackScore;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained) return trk.chi2pid[trk.bestplane].chi2_muon;
        else return -999.;

      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained) return trk.chi2pid[trk.bestplane].chi2_proton;
        else return -999.;

      }
      else{
        return -999.;
      }
    });

  } // end namespace Aux


} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
