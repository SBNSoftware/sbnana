#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  const Cut HasTwoPrimaryTracks([](const caf::SRSliceProxy* slc) {
    return ICARUSNumuXsec::PrimaryTrackIndices(slc).size()>=2;
  });
  const Cut HasMuonTrack([](const caf::SRSliceProxy* slc) {
    return MuonTrackIndex(slc)>=0.;
  });

  const Cut HasProtonTrack([](const caf::SRSliceProxy* slc) {
    return ProtonTrackIndex(slc)>=0;
  });

  const Cut MuonTrackContained([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      return isContained;
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackExiting([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      return !isContained;
    }
    else{
      return false;
    }
  });

  const Cut MuonTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      if(intID==-1) return false;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      return ( abs(trk_truth_pdg)==13 && isTrkTruthContained );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthExitingNuMuon([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      if(intID==-1) return false;
      
      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      return ( abs(trk_truth_pdg)==13 && !isTrkTruthContained );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthCosmicMuon([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      return (intID==-1);
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthStoppingProton([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isProtonStopped = (endProc==2);
      return ( trk_truth_pdg==2212 && isTrkTruthContained && isProtonStopped );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthInelProton([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isProtonInel = (endProc==7);
      return ( trk_truth_pdg==2212 && isTrkTruthContained && isProtonInel );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthOtherProton([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isProtonStopped = (endProc==2);
      bool isProtonInel = (endProc==7);

      return ( trk_truth_pdg==2212 && 
               ( !isTrkTruthContained ||
                 (isTrkTruthContained && !(isProtonStopped||isProtonInel))
               )
             );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthStoppingChargedPion([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      //int intID = trk_truth.interaction_id;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      return ( abs(trk_truth_pdg)==211 && isTrkTruthContained );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthExitingChargedPion([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      //int intID = trk_truth.interaction_id;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      return ( abs(trk_truth_pdg)==211 && !isTrkTruthContained );
    }
    else{
      return false;
    }
  });

  const Cut MuonTrackTruthOther([](const caf::SRSliceProxy* slc) {
    bool IsCategorized = MuonTrackTruthContainedNuMuon(slc);
    IsCategorized |= MuonTrackTruthExitingNuMuon(slc);
    IsCategorized |= MuonTrackTruthCosmicMuon(slc);
    IsCategorized |= MuonTrackTruthStoppingProton(slc);
    IsCategorized |= MuonTrackTruthInelProton(slc);
    IsCategorized |= MuonTrackTruthOtherProton(slc);
    IsCategorized |= MuonTrackTruthStoppingChargedPion(slc);
    IsCategorized |= MuonTrackTruthExitingChargedPion(slc);
    return !IsCategorized;
  });


  const Cut ProtonPCut([](const caf::SRSliceProxy* slc) {
    double protonP = ProtonTrackP(slc);
    if(protonP<0){
      return false;
    }
    else{
      return (protonP>0.6);
    }
  });

  // pion tagging
  const Cut NoPionCand([](const caf::SRSliceProxy* slc) {

    int muonTrackIndex = MuonTrackIndex(slc);
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(muonTrackIndex<0 || protonTrackIndex<0) return false;
/*
    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
      const auto& pfp = slc->reco.pfp.at(i_pfp);
      const auto& trk = pfp.trk;

      if(i_pfp==(unsigned int)muonTrackIndex) continue;
      if(i_pfp==(unsigned int)protonTrackIndex) continue;

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
*/
    return false;

  });

  namespace Aux{

    const Cut HasRelaxedMuonTrack([](const caf::SRSliceProxy* slc) {
      return RelaxedMuonTrackIndex(slc)>=0;
    });
    const Cut HasRelaxedProtonTrack([](const caf::SRSliceProxy* slc) {
      return RelaxedProtonTrackIndex(slc)>=0;
    });

    const Cut RelaxedMuonTrackContained([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(isnan(trk.end.x) || isnan(trk.end.y) || isnan(trk.end.z)) return false;
        bool isTrkContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        return isTrkContained;
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackContained([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(isnan(trk.end.x) || isnan(trk.end.y) || isnan(trk.end.z)) return false;
        bool isTrkContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        return isTrkContained;
      }
      else{
        return false;
      }
    });

    const Cut RelaxedMuonTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthExitingNuMuon([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;
        
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && !isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthCosmicMuon([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        return (intID==-1);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthStoppingProton([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isProtonStopped = (endProc==2);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isProtonStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthInelProton([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isProtonInel = (endProc==7);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isProtonInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthOtherProton([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isProtonStopped = (endProc==2);
        bool isProtonInel = (endProc==7);

        return ( trk_truth_pdg==2212 && 
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isProtonStopped||isProtonInel))
                 )
               );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthStoppingChargedPion([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        //int intID = trk_truth.interaction_id;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthExitingChargedPion([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        //int intID = trk_truth.interaction_id;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==211 && !isTrkTruthContained );
      }
      else{
        return false;
      }
    });

    const Cut RelaxedMuonTrackTruthOther([](const caf::SRSliceProxy* slc) {
      bool IsCategorized = RelaxedMuonTrackTruthContainedNuMuon(slc);
      IsCategorized |= RelaxedMuonTrackTruthExitingNuMuon(slc);
      IsCategorized |= RelaxedMuonTrackTruthCosmicMuon(slc);
      IsCategorized |= RelaxedMuonTrackTruthStoppingProton(slc);
      IsCategorized |= RelaxedMuonTrackTruthInelProton(slc);
      IsCategorized |= RelaxedMuonTrackTruthOtherProton(slc);
      IsCategorized |= RelaxedMuonTrackTruthStoppingChargedPion(slc);
      IsCategorized |= RelaxedMuonTrackTruthExitingChargedPion(slc);
      return !IsCategorized;
    });


    const Cut RelaxedProtonTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthExitingNuMuon([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;
        
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && !isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthCosmicMuon([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        return (intID==-1);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthStoppingProton([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isProtonStopped = (endProc==2);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isProtonStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthInelProton([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;
        
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isProtonInel = (endProc==7);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isProtonInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthOtherProton([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isProtonStopped = (endProc==2);
        bool isProtonInel = (endProc==7);

        return ( trk_truth_pdg==2212 &&
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isProtonStopped||isProtonInel))
                 )
               );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthStoppingChargedPion([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthExitingChargedPion([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
      
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==211 && !isTrkTruthContained );
      }
      else{
        return false;
      } 
    });
    const Cut RelaxedProtonTrackTruthOther([](const caf::SRSliceProxy* slc) {
      bool IsCategorized = RelaxedProtonTrackTruthContainedNuMuon(slc);
      IsCategorized |= RelaxedProtonTrackTruthExitingNuMuon(slc);
      IsCategorized |= RelaxedProtonTrackTruthCosmicMuon(slc);
      IsCategorized |= RelaxedProtonTrackTruthStoppingProton(slc);
      IsCategorized |= RelaxedProtonTrackTruthInelProton(slc);
      IsCategorized |= RelaxedProtonTrackTruthOtherProton(slc);
      IsCategorized |= RelaxedProtonTrackTruthStoppingChargedPion(slc);
      IsCategorized |= RelaxedProtonTrackTruthExitingChargedPion(slc);
      return !IsCategorized;
    });

    const Cut RelaxedProtonTrackIsochronous([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return ( fabs(trk.dir.x)<0.1 );
       
      }
      else{
        return false;
      }
    });

    const Cut RelaxedProtonTrackHasChi2MuonExcess([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return false;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        return (Chi2Muon>40. && Chi2Muon<50.);
      }
      else{
        return false;
      }
    });

    const Cut RelaxedProtonTrackHasChi2MuonDeficit([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return false;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        return (Chi2Muon>50.);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackShort([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return (trk.len<10.);
      }
      else{
        return false;
      }
    });

  } // end namespace Aux


} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
