#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana {


// Mostly just copied from NumuVarsIcarus202106.cxx
const Var kIcarus202401MuonIdx([](const caf::SRSliceProxy* slc) -> int {
      // The (dis)qualification of a slice is based upon the track level information.
      float Longest(0);
      int PTrackInd(-1);
      for (std::size_t i(0); i < slc->reco.npfp; ++i) {
        auto const& pfp = slc->reco.pfp.at(i);
        if (pfp.trackScore < 0.5) { continue; }
        auto const& trk = pfp.trk;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

        int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;

        float Chi2Proton = trk.chi2pid[plane].chi2_proton;
        float Chi2Muon = trk.chi2pid[plane].chi2_muon;
  
        const bool Contained = ( !isnan(trk.end.x) &&
            (( trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5 ) ||
             ( trk.end.x > 61.94 + 5 && trk.end.x < +358.49 - 5 )) &&
            !isnan(trk.end.y) &&
            ( trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5 ) &&
            !isnan(trk.end.z) &&
            ( trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5 ) );
        const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest ) {
          Longest = trk.len;
          PTrackInd = i;
        }
      }
      return PTrackInd;
});

static bool Icarus202401_proton_cut(const caf::SRTrackProxy& trk)
{
  int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
  return trk.chi2pid[plane].chi2_proton < 100 && trk.chi2pid[plane].chi2_proton > 0;
}

const Var kIcarus202401NumPions([](const caf::SRSliceProxy* slc)
{
  int count = 0;
  auto idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    if (pfp.trackScore < 0.5) { continue; }
    auto const& trk = pfp.trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
    if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
    if(pfp.id != muID && !Icarus202401_proton_cut(trk) && AtSlice && std::hypot(pfp.trk.rangeP.p_pion, 0.140f) - 0.14 > 0.025)
      ++count;
  }
  return count;
});

const Var kIcarus202401NumProtons([](const caf::SRSliceProxy* slc)
{
  int count = 0;
  auto idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if (idx >= 0) muID = slc->reco.pfp.at(idx).id;
  for(auto& pfp: slc->reco.pfp) {
    if (pfp.trackScore < 0.5) { continue; }
    if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
    auto const& trk = pfp.trk;
    const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                   slc->vertex.y - trk.start.y,
                                   slc->vertex.z - trk.start.z);
    const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
    if(pfp.id != muID && Icarus202401_proton_cut(trk) && AtSlice && std::hypot(pfp.trk.rangeP.p_proton, 0.938f) - 0.938 > 0.025)
      ++count;
  }
  return count;
});

const Var kIcarus202401NumShowers([](const caf::SRSliceProxy* slc){
  int count = 0;
  for(const auto& pfp: slc->reco.pfp) {
    const float Atslc = std::hypot(slc->vertex.x - pfp.shw.start.x,
                                   slc->vertex.y - pfp.shw.start.y,
                                   slc->vertex.z - pfp.shw.start.z);
    const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
    if(pfp.trackScore > 0 && pfp.trackScore < 0.5 && AtSlice && pfp.shw.plane[2].energy > 0.025) count++;
  }
  return count;
});

const Var kIcarus202401RecoMuonE([](const caf::SRSliceProxy* slc){
  return std::hypot(kIcarus202401RecoMuonP(slc), 0.106f);
});

const Var kIcarus202401RecoProtonKE([](const caf::SRSliceProxy* slc){
  double E = 0;
  for(const auto P: kIcarus202401RecoProtonP(slc))
    E += std::hypot(P, 0.938f) - 0.938;
  return E;
});

const Var kIcarus202401RecoPionE([](const caf::SRSliceProxy* slc){
  double E = 0;
  for(const auto P: kIcarus202401RecoPionP(slc))
    E += std::hypot(P, 0.140f);
  return E;
});

const Var kIcarus202401RecoNuE = kIcarus202401RecoMuonE + kIcarus202401RecoProtonKE + kIcarus202401RecoPionE;

static bool Icarus202401_contained(const caf::SRTrackProxy& trk)
{
  return ((trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5)
             || (trk.end.x > 61.94 + 5 && trk.end.x < 358.49 - 5))
           && trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5
           && trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5;
}

const Var kIcarus202401RecoMuonP([](const caf::SRSliceProxy *slc){
  double ret = 0;
  if(kIcarus202401MuonIdx(slc) < 0) return ret;
  const auto &muTrk = slc->reco.pfp.at(kIcarus202401MuonIdx(slc)).trk;
  ret = Icarus202401_contained(muTrk) ? muTrk.rangeP.p_muon : 
               muTrk.mcsP.is_bwd_muon ? muTrk.mcsP.bwdP_muon : 
                                        muTrk.mcsP.fwdP_muon;
  return ret;
});

const MultiVar kIcarus202401RecoProtonP([](const caf::SRSliceProxy* slc){
  std::vector<double> Ps;
  int mu_idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if(mu_idx >= 0) muID = slc->reco.pfp[mu_idx].id;
  for(const auto& pfp: slc->reco.pfp) {
    if(pfp.id == muID) continue;
    if(pfp.trackScore < 0.5) continue;
    if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
    if(Icarus202401_proton_cut(pfp.trk) && pfp.parent_is_primary)
      Ps.push_back(pfp.trk.rangeP.p_proton);
  }
  return Ps;
});

const MultiVar kIcarus202401RecoPionP([](const caf::SRSliceProxy* slc){
  std::vector<double> Ps;
  int mu_idx = kIcarus202401MuonIdx(slc);
  int muID = -1;
  if(mu_idx >= 0) muID = slc->reco.pfp[mu_idx].id;
  for(const auto& pfp: slc->reco.pfp) {
    if(pfp.id == muID) continue;
    if(pfp.trackScore < 0.5) continue;
    if(pfp.trk.chi2pid[2].chi2_proton == 0 && pfp.trk.chi2pid[1].chi2_proton == 0) continue;
    if(!Icarus202401_proton_cut(pfp.trk) && pfp.parent_is_primary)
      Ps.push_back(pfp.trk.rangeP.p_proton);
  }
  return Ps;
});
}
