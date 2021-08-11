#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana {
// Get the Largest shower in the slice
const Var kLargestRecoShowerIdx(
    [](const caf::SRSliceProxy* slc) -> int {
      int bestIdx(-1);
      double maxEnergy(-1);

      for (unsigned int i = 0; i < slc->reco.nshw; i++) {
        auto const& shw = slc->reco.shw[i];
        if (!shw.pfp.parent_is_primary)
          continue;
        if (shw.bestplane_energy > maxEnergy) {
          bestIdx = i;
          maxEnergy = shw.bestplane_energy;
        }
      }
      return bestIdx;
    });

// Currently assumes shw 0 is the primary
const Var kRecoShower_BestEnergy(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].bestplane_energy;
    });

const Var kRecoShower_TruePdg(
    [](const caf::SRSliceProxy* slc) -> int {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].truth.p.pdg;
    });

// Currently assumes shw 0 is the primary
const Var kRecoShower_BestdEdx(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      const double dEdx = slc->reco.shw[largestShwIdx].bestplane_dEdx;
      if(dEdx < 0) return -5;
      return dEdx;
    });

const Var kRecoShower_ConversionGap(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].conversion_gap;
    });

const Var kRecoShower_Density(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].density;
    });

const Var kRecoShower_Energy(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx != -1) return -5;
      return slc->reco.shw[largestShwIdx].energy[1]; // so far taking whatever plane 1 is and first shw
    });

const Var kRecoShower_Length(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].len;
    });

const Var kRecoShower_OpenAngle(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return 180 * slc->reco.shw[largestShwIdx].open_angle / M_PI;
    });

const Var kRecoShower_StartX(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999;
      return slc->reco.shw[largestShwIdx].start.x;
    });

const Var kRecoShower_StartY(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx != -1) return -9999;
      return slc->reco.shw[largestShwIdx].start.y;
    });

const Var kRecoShower_StartZ(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999;
      return slc->reco.shw[largestShwIdx].start.z;
    });

const Var kRecoShower_EndX(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999;
      const caf::SRShowerProxy& shw = slc->reco.shw[largestShwIdx];
      return shw.start.x + shw.dir.x * shw.len;
    });

const Var kRecoShower_EndY(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999;
      const caf::SRShowerProxy& shw = slc->reco.shw[largestShwIdx];
      return shw.start.y + shw.dir.y * shw.len;
    });

const Var kRecoShower_EndZ(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999;
      const caf::SRShowerProxy& shw = slc->reco.shw[largestShwIdx];
      return shw.start.z + shw.dir.z * shw.len;
    });

const Var kRecoShower_densityGradient(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].selVars.densityGradient;
    });

const Var kRecoShower_densityGradientPower(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].selVars.densityGradientPower;
    });

const Var kRecoShower_trackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].selVars.trackLength;
    });

const Var kRecoShower_trackWidth(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5;
      return slc->reco.shw[largestShwIdx].selVars.trackWidth;
    });

const Var kRecoShowers_EnergyCut(
    [](const caf::SRSliceProxy* slc) -> double {
      unsigned int counter(0);
      for (auto const& shw : slc->reco.shw) {
        if (shw.bestplane_energy > 200)
          ++counter;
      }
      return counter;
    });

const Var kLongestTrackIdx(
    [](const caf::SRSliceProxy* slc) -> int {
      int bestIdx(-1);
      double maxLength(-1);

      for (unsigned int i = 0; i < slc->reco.ntrk; i++) {
        auto const& trk = slc->reco.trk[i];
        if (!trk.pfp.parent_is_primary)
          continue;

        if (trk.len > maxLength) {
          bestIdx = i;
          maxLength = trk.len;
        }
      }
      return bestIdx;
    });

const Var kLongestTrackTruePdg(
    [](const caf::SRSliceProxy* slc) -> int {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      if(longestTrackIdx != -1) return -5;
      return slc->reco.trk[longestTrackIdx].truth.p.pdg;
    });

const Var kLongestTrackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      if(longestTrackIdx != -1) return -5;
      return slc->reco.trk[longestTrackIdx].len;
    });

const Var kLongestTrackChi2Muon(
    [](const caf::SRSliceProxy* slc) -> double {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      float chi2(-5.f);
      if (longestTrackIdx != -1) {
        const unsigned int bestplane(slc->reco.trk[longestTrackIdx].bestplane);

        chi2 = slc->reco.trk[longestTrackIdx].chi2pid[bestplane].chi2_muon;
      }
      return chi2;
    });

const Var kLongestTrackChi2Pion(
    [](const caf::SRSliceProxy* slc) -> double {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      float chi2(-5.f);
      if (longestTrackIdx != -1) {
        const unsigned int bestplane(slc->reco.trk[longestTrackIdx].bestplane);

        chi2 = slc->reco.trk[longestTrackIdx].chi2pid[bestplane].chi2_pion;
      }
      return chi2;
    });

const Var kLongestTrackChi2Kaon(
    [](const caf::SRSliceProxy* slc) -> double {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      float chi2(-5.f);
      if (longestTrackIdx != -1) {
        const unsigned int bestplane(slc->reco.trk[longestTrackIdx].bestplane);

        chi2 = slc->reco.trk[longestTrackIdx].chi2pid[bestplane].chi2_kaon;
      }
      return chi2;
    });

const Var kLongestTrackChi2Proton(
    [](const caf::SRSliceProxy* slc) -> double {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      float chi2(-5.f);
      if (longestTrackIdx != -1) {
        const unsigned int bestplane(slc->reco.trk[longestTrackIdx].bestplane);

        chi2 = slc->reco.trk[longestTrackIdx].chi2pid[bestplane].chi2_proton;
      }
      return chi2;
    });

const Var kMuonTrackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      double length = -5.f;
      const int longestTrackIdx(kLongestTrackIdx(slc));
      if(longestTrackIdx != -1 && (kLongestTrackChi2Muon(slc) < 30.f && kLongestTrackChi2Proton(slc) > 60.f)) {
        length = slc->reco.trk[longestTrackIdx].len;
      }
      return length;
    });

const Var kLongestTrackDazzlePID(
    [](const caf::SRSliceProxy* slc) -> int {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      return longestTrackIdx != -1 ? (int)slc->reco.trk[longestTrackIdx].dazzle.pdg : -5;
    });

const Var kLongestTrackDazzleMuonScore(
    [](const caf::SRSliceProxy* slc) -> float {
      const int longestTrackIdx(kLongestTrackIdx(slc));
      return longestTrackIdx != -1 ? (float)slc->reco.trk[longestTrackIdx].dazzle.muonScore : -5.f;
    });

const Var kRecoShowerRazzlePID(
    [](const caf::SRSliceProxy* slc) -> int {
      const int showerIdx(kLargestRecoShowerIdx(slc));
      return showerIdx != -1 ? (int)slc->reco.shw[showerIdx].razzle.pdg : -5;
    });

const Var kRecoShowerRazzleElectronScore(
    [](const caf::SRSliceProxy* slc) -> float {
      const int showerIdx(kLargestRecoShowerIdx(slc));
      return showerIdx != -1 ? (float)slc->reco.shw[showerIdx].razzle.electronScore : -5.f;
    });

}
