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

        // pid from collection
        const float Chi2Proton = trk.chi2pid[2].chi2_proton;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

        const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
        const bool MaybeMuonContained = ( Contained && trk.calo[2].nhit!=0 && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
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
  const Var MuonTrackP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        if(isnan(trk.rangeP.p_muon)) return -999.;
        return trk.rangeP.p_muon;
      }
      else{
        if(isnan(trk.mcsP.fwdP_muon)) return -999.;
        return trk.mcsP.fwdP_muon;
      }
    }
    else{
      return -999.;
    }

  });
  const Var MuonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
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
  // truth match
  const Var MuonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      if( isnan(trk.truth.p.length) ) return -999.;
      return trk.truth.p.length;
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      if( isnan(trk.truth.p.genE) ) return -999.;
      return trk.truth.p.genE - M_MUON;
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& truth_p = trk.truth.p;
      if( isnan(truth_p.genp.x) ) return -999.;
      TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
      double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });

  // proton
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

        if ( AtSlice && Contained && trk.calo[2].nhit!=0 && Chi2Proton <= 100 && Chi2Muon >= 30 && angle >= -0.9 && trk.len > Longest ) {
          Longest = trk.len;
          PTrackInd = trkIdx;
        }

      }

      return PTrackInd;

    }
  });

  const Var ProtonTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      return trk.len;
    }
    else{
      return -999.;
    }

  });
  const Var ProtonTrackP([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        if(isnan(trk.rangeP.p_proton)) return -999.;
        return trk.rangeP.p_proton;
      }
      else return -999.;
    }
    else{
      return -999.;
    }

  });
  const Var ProtonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
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
      static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      
      TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;
      
      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  // truth match
  const Var ProtonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      if( isnan(trk.truth.p.length) ) return -999.;
      return trk.truth.p.length;
    }
    else{
      return -999.;
    }
  }); 
  const Var ProtonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      if( isnan(trk.truth.p.genE) ) return -999.;
      return trk.truth.p.genE - M_PROTON;
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& truth_p = trk.truth.p;
      if( isnan(truth_p.genp.x) ) return -999.;
      TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
      double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });

  // muon+proton
  const Var MuonProtonCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(muonTrackIndex>=0 && protonTrackIndex>=0){
      const auto& trk_mu = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_pro = slc->reco.pfp.at(protonTrackIndex).trk;

      TVector3 v3_trk_mu(trk_mu.dir.x, trk_mu.dir.y, trk_mu.dir.z);
      TVector3 v3_trk_pro(trk_pro.dir.x, trk_pro.dir.y, trk_pro.dir.z);

      double angle = v3_trk_mu.Angle(v3_trk_pro);
      return TMath::Cos(angle);

    }
    else{
      return -999.;
    }
  });

  namespace Aux{

    const Var RelaxedMuonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
      if(slc->reco.pfp.size()==0){
        return -999.;
      }
      else{
        int muonIdx = -1;
        double muonMaxLength = -1;
        for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
          const auto& pfp = slc->reco.pfp.at(i_pfp);
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

          if ( AtSlice && trk.len > muonMaxLength ) {
            muonIdx = i_pfp;
            muonMaxLength = trk.len;
          }
        }
        return muonIdx;

      }
    });

    const Var RelaxedMuonTrackLength([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.len;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackP([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained){
          if(isnan(trk.rangeP.p_muon)) return -999.;
          return trk.rangeP.p_muon;
        }
        else{
          if(isnan(trk.mcsP.fwdP_muon)) return -999.;
          return trk.mcsP.fwdP_muon;
        }

      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
        double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackDirX([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.dir.x;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackDirY([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.dir.y;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackDirZ([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.dir.z;
      }
      else{
        return -999.;
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

        return trk.chi2pid[trk.bestplane].chi2_muon;
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

        return trk.chi2pid[trk.bestplane].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;
        
        return trk.chi2pid[2].chi2_muon;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;
        
        return trk.chi2pid[2].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackCustomChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 3);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackCustomChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 0);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const MultiVar RelaxedMuonTrackCollectionRR([](const caf::SRSliceProxy* slc) -> vector<double> {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      vector<double> rets;
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.rr );
        }
      }
      return rets;
    });
    const MultiVar RelaxedMuonTrackCollectiondEdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      vector<double> rets;
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dedx );
        }
      }
      return rets;
    });
    const MultiVar RelaxedMuonTrackCollectiondQdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      vector<double> rets;
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dqdx );
        }
      }
      return rets;
    });
    // truth match
    const Var RelaxedMuonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if( isnan(trk.truth.p.length) ) return -999.;
        return trk.truth.p.length;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedMuonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if( isnan(trk.truth.p.genE) ) return -999.;
        double P2 = trk.truth.p.genE*trk.truth.p.genE - M_MUON*M_MUON;
        if(P2>0) return std::sqrt(P2);
        else return -999.;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthPResFrac([](const caf::SRSliceProxy* slc) -> double {
      double Reco = RelaxedMuonTrackP(slc);
      double True = RelaxedMuonTrackTruthP(slc);
      if(Reco>0. && True>0.){
        return (Reco-True)/True;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthOneOverP([](const caf::SRSliceProxy* slc) -> double {
      double P = RelaxedMuonTrackTruthP(slc);
      if(P>0.) return 1./P;
      else return -999.;
    });
    const Var RelaxedMuonTrackTruthOneOverPResFrac([](const caf::SRSliceProxy* slc) -> double {
      double RecoP = RelaxedMuonTrackP(slc);
      double True = RelaxedMuonTrackTruthOneOverP(slc);
      if(RecoP>0. && True>0.){
        double Reco = 1./RecoP;
        return (Reco-True)/True;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedMuonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& truth_p = trk.truth.p;
        if( isnan(truth_p.genp.x) ) return -999.;
        TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
        double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthStartProcess([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        return trk.truth.p.start_process;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthEndProcess([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        return trk.truth.p.end_process;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedProtonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
      if(slc->reco.pfp.size()==0){
        return -999.;
      }
      else{

        // Requiring good muon
        int muonTrackIndex = MuonTrackIndex(slc);
        int protonIdx = -1;

        const int chargedpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
        const bool HasPionCand = (chargedpionTrackIndex>=0);

        if(muonTrackIndex>=0 && !HasPionCand){

          double protonMaxLength = -1.;

          // find the longest track

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

            const unsigned int idxPrim = (unsigned int)muonTrackIndex;
            const TVector3 muDir( slc->reco.pfp[idxPrim].trk.dir.x, slc->reco.pfp[idxPrim].trk.dir.y, slc->reco.pfp[idxPrim].trk.dir.z );
            const TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
            const float angle = TMath::Cos(muDir.Angle(pDir));

            if( AtSlice && angle>=-0.9 && trk.len > protonMaxLength ){
              protonMaxLength = trk.len;
              protonIdx = i_pfp;
            }


          }

        }

        if(HasPionCand) return -1;
        return protonIdx;

      }
    });
    const Var RelaxedProtonTrackLength([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.len;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackP([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained){
          if(isnan(trk.rangeP.p_proton)) return -999.;
          return trk.rangeP.p_proton;
        }
        else{
          if(isnan(trk.mcsP.fwdP_proton)) return -999.;
          return trk.mcsP.fwdP_proton;
        }

      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
        double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackDirX([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.x;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackDirY([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.y;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackDirZ([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.z;
      }
      else{
        return -999.;
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

        return trk.chi2pid[trk.bestplane].chi2_muon;
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

        return trk.chi2pid[trk.bestplane].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;
        return trk.chi2pid[2].chi2_muon;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;

        return trk.chi2pid[2].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackCustomChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 3);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackCustomChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 0);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackNHitCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return trk.calo[2].nhit;
      }
      else{
        return -999.;
      }
    });
    const MultiVar RelaxedProtonTrackCollectionRR([](const caf::SRSliceProxy* slc) -> vector<double> {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      vector<double> rets;
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.rr );
        }
      }
      return rets;
    });
    const MultiVar RelaxedProtonTrackCollectiondEdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      vector<double> rets;
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dedx );
        }
      }
      return rets;
    });
    const MultiVar RelaxedProtonTrackCollectiondQdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      vector<double> rets;
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dqdx );
        }
      }
      return rets;
    });
    // truth match
    const Var RelaxedProtonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if( isnan(trk.truth.p.length) ) return -999.;
        return trk.truth.p.length;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if( isnan(trk.truth.p.genE) ) return -999.;
        double P2 = trk.truth.p.genE*trk.truth.p.genE - M_PROTON*M_PROTON;
        if(P2>0) return std::sqrt(P2);
        else return -999.;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthPResFrac([](const caf::SRSliceProxy* slc) -> double {
      double Reco = RelaxedProtonTrackP(slc);
      double True = RelaxedProtonTrackTruthP(slc);
      if(Reco>0. && True>0.){
        return (Reco-True)/True;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& truth_p = trk.truth.p;
        if( isnan(truth_p.genp.x) ) return -999.;
        TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
        double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthStartProcess([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return trk.truth.p.start_process;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthEndProcess([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return trk.truth.p.end_process;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedChargedPionTrackIndex([](const caf::SRSliceProxy* slc) -> double {
      vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
      if(primTrackIndices.size()==0){
        return -999.;
      }
      else{

        // Requiring good muon
        int muonTrackIndex = MuonTrackIndex(slc);
        int PTrackInd(-1);
        if(muonTrackIndex>=0){

          float Longest(0);
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

            const float Chi2Proton = trk.chi2pid[2].chi2_proton;
            const float Chi2Muon = trk.chi2pid[2].chi2_muon;
            const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
            const bool MaybeTrackExiting = ( !Contained && trk.len > 100);
            const bool MaybeTrackContained = ( Contained && trk.calo[2].nhit!=0 && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 20 );

            // check pion cand other than the good muon we found
            if ( AtSlice && ( MaybeTrackExiting || MaybeTrackContained ) && trk.len > Longest ){
              Longest = trk.len;
              PTrackInd = i_pfp;
            }
          } // END pfp loop

        } // END if muon track exist

        return PTrackInd;

      } // END if slice has primary tracks
    });


  } // end namespace Aux


} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
