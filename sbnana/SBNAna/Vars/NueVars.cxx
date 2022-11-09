#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana {
// Get the Largest shower in the slice
const Var kLargestRecoShowerIdx(
    [](const caf::SRSliceProxy* slc) -> int {
      int bestIdx(-1);
      double maxEnergy(-1);

      for (unsigned int i = 0; i < slc->reco.npfp; i++) {
        auto const& shw = slc->reco.pfp[i].shw;
        if (!slc->reco.pfp[i].parent_is_primary || slc->reco.pfp[i].trackScore > 0.5)
          continue;
        if (shw.bestplane_energy > maxEnergy) {
          bestIdx = i;
          maxEnergy = shw.bestplane_energy;
        }
      }
      return bestIdx;
    });

/// Pointer to largest reconstructed shower, or null pointer if none exists
const caf::SRShowerProxy* LargestRecoShower(const caf::SRSliceProxy* slc)
{
  const int largestShwIdx = kLargestRecoShowerIdx(slc);
  if(largestShwIdx == -1) return 0;
  return &slc->reco.pfp[largestShwIdx].shw;
}

// Currently assumes shw 0 is the primary
const Var kRecoShower_BestEnergy(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->bestplane_energy) : -5;
    });

const Var kRecoShower_TruePdg(
    [](const caf::SRSliceProxy* slc) -> int {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? int(shw->truth.p.pdg) : -5;
    });

// Currently assumes shw 0 is the primary
const Var kRecoShower_BestdEdx(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      if(!shw || shw->bestplane_dEdx < 0) return -5.;
      return shw->bestplane_dEdx;
    });

const Var kRecoShower_ConversionGap(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->conversion_gap) : -5.;
    });

const Var kRecoShower_Density(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->density) : -5.;
    });

const Var kRecoShower_Energy(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->plane[1].energy) : -5.; // so far taking whatever plane 1 is and first shw
    });

const Var kRecoShower_Length(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->len) : -5.;
    });

const Var kRecoShower_OpenAngle(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? 180. * shw->open_angle / M_PI : -5.;
    });

const Var kRecoShower_StartX(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->start.x) : -9999.;
    });

const Var kRecoShower_StartY(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->start.y) : -9999.;
    });

const Var kRecoShower_StartZ(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->start.z) : -9999.;
    });

const Var kRecoShower_EndX(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->start.x + shw->dir.x * shw->len) : -9999.;
    });

const Var kRecoShower_EndY(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->start.y + shw->dir.y * shw->len) : -9999.;
    });

const Var kRecoShower_EndZ(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->start.z + shw->dir.z * shw->len) : -9999.;
    });

const Var kRecoShower_densityGradient(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->selVars.densityGradient) : -5.;
    });

const Var kRecoShower_densityGradientPower(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->selVars.densityGradientPower) : -5.;
    });

const Var kRecoShower_trackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->selVars.trackLength) : -5.;
    });

const Var kRecoShower_trackWidth(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? double(shw->selVars.trackWidth) : -5.;
    });

const Var kRecoShowers_EnergyCut(
    [](const caf::SRSliceProxy* slc) -> unsigned {
      unsigned int counter(0);
      for (auto const& pfp : slc->reco.pfp) {
        if (pfp.trackScore > 0.5) { continue; }
	auto const& shw = pfp.shw;
        if (shw.bestplane_energy > 0.2f)
          ++counter;
      }
      return counter;
    });

const Var kLongestTrackIdx(
    [](const caf::SRSliceProxy* slc) -> int {
      int bestIdx(-1);
      double maxLength(-1);

      for (unsigned int i = 0; i < slc->reco.npfp; i++) {
        auto const& trk = slc->reco.pfp[i].trk;
        if (!slc->reco.pfp[i].parent_is_primary || slc->reco.pfp[i].trackScore < 0.5)
          continue;

        if (trk.len > maxLength) {
          bestIdx = i;
          maxLength = trk.len;
        }
      }
      return bestIdx;
    });

/// Pointer to longest reconstructed shower, or null pointer if none exists
const caf::SRTrackProxy* LongestRecoTrack(const caf::SRSliceProxy* slc)
{
  const int longestTrackIdx = kLongestTrackIdx(slc);
  if(longestTrackIdx == -1) return 0;
  return &slc->reco.pfp[longestTrackIdx].trk;
}

const Var kLongestTrackTruePdg(
    [](const caf::SRSliceProxy* slc) -> int {
      const caf::SRTrackProxy* trk = LongestRecoTrack(slc);
      return trk ? int(trk->truth.p.pdg) : -5;
    });

const Var kLongestTrackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRTrackProxy* trk = LongestRecoTrack(slc);
      return trk ? double(trk->len) : -5.;
    });


const caf::SRTrkChi2PIDProxy* BestPlaneChi2PID(const caf::SRTrackProxy* trk)
{
  if(trk->bestplane == -1) return 0;
  return &trk->chi2pid[trk->bestplane];
}

const caf::SRTrkChi2PIDProxy* LongestTrackBestPlaneChi2PID(const caf::SRSliceProxy* slc)
{
  const caf::SRTrackProxy* trk = LongestRecoTrack(slc);
  return trk ? BestPlaneChi2PID(trk) : 0;
}

const Var kLongestTrackChi2Muon(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRTrkChi2PIDProxy* chi2 = LongestTrackBestPlaneChi2PID(slc);
      return chi2 ? double(chi2->chi2_muon) : -5.;
    });

const Var kLongestTrackChi2Pion(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRTrkChi2PIDProxy* chi2 = LongestTrackBestPlaneChi2PID(slc);
      return chi2 ? double(chi2->chi2_pion) : -5.;
    });

const Var kLongestTrackChi2Kaon(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRTrkChi2PIDProxy* chi2 = LongestTrackBestPlaneChi2PID(slc);
      return chi2 ? double(chi2->chi2_kaon) : -5.;
    });

const Var kLongestTrackChi2Proton(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRTrkChi2PIDProxy* chi2 = LongestTrackBestPlaneChi2PID(slc);
      return chi2 ? double(chi2->chi2_proton) : -5.;
    });

const Var kMuonTrackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRTrackProxy* trk = LongestRecoTrack(slc);
      if(!trk) return -5;
      if(trk && (kLongestTrackChi2Muon(slc) < 30.f && kLongestTrackChi2Proton(slc) > 60.f)) {
        return trk->len;
      }
      return -5;
    });

const Var kLongestTrackDazzlePID(
    [](const caf::SRSliceProxy* slc) -> int {
      const caf::SRTrackProxy* trk = LongestRecoTrack(slc);
      return trk ? (int)trk->dazzle.pdg : -5;
    });

const Var kLongestTrackDazzleMuonScore(
    [](const caf::SRSliceProxy* slc) -> float {
      const caf::SRTrackProxy* trk = LongestRecoTrack(slc);
      return trk ? (float)trk->dazzle.muonScore : -5.f;
    });

const Var kRecoShowerRazzlePID(
    [](const caf::SRSliceProxy* slc) -> int {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? (int)shw->razzle.pdg : -5;
    });

const Var kRecoShowerRazzleElectronScore(
    [](const caf::SRSliceProxy* slc) -> float {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      return shw ? (float)shw->razzle.electronScore : -5.f;
    });

}
