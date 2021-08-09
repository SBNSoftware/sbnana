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

const caf::SRShowerProxy* LargestRecoShower(const caf::SRSliceProxy* slc)
{
  const int largestShwIdx(kLargestRecoShowerIdx(slc));
  if(largestShwIdx == -1) return 0;
  return &slc->reco.shw[largestShwIdx];
}

// Currently assumes shw 0 is the primary
const Var kRecoShower_BestEnergy(
    [](const caf::SRSliceProxy* slc) -> double {
      const caf::SRShowerProxy* shw = LargestRecoShower(slc);
      if(!shw) return -5.;
      return shw->bestplane_energy;
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
      if(largestShwIdx == -1) return -5.;
      const double dedx = slc->reco.shw[largestShwIdx].bestplane_dEdx;
      if(dedx > 0) return dedx; else return -5.;
    });

const Var kRecoShower_ConversionGap(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5.;
      return slc->reco.shw[largestShwIdx].conversion_gap;
    });

const Var kRecoShower_Density(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5.;
      return slc->reco.shw[largestShwIdx].density;
    });

const Var kRecoShower_Energy(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5.;
      return slc->reco.shw[largestShwIdx].energy[1]; // so far taking whatever plane 1 is and first shw
    });

const Var kRecoShower_Length(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5.;
      return slc->reco.shw[largestShwIdx].len;
    });

const Var kRecoShower_OpenAngle(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -5.;
      return 180 * slc->reco.shw[largestShwIdx].open_angle / M_PI;
    });

const Var kRecoShower_StartX(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999.;
      return slc->reco.shw[largestShwIdx].start.x;
    });

const Var kRecoShower_StartY(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999.;
      return slc->reco.shw[largestShwIdx].start.y;
    });

const Var kRecoShower_StartZ(
    [](const caf::SRSliceProxy* slc) -> double {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if(largestShwIdx == -1) return -9999.;
      return slc->reco.shw[largestShwIdx].start.z;
    });

const Var kRecoShower_EndX(
    [](const caf::SRSliceProxy* slc) -> double {
      double end = -9999.0;
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if (largestShwIdx != -1) {
        double start = slc->reco.shw[largestShwIdx].start.x;
        double dir = slc->reco.shw[largestShwIdx].dir.x;
        double len = slc->reco.shw[largestShwIdx].len;
        end = start + (dir * len);
      }
      return end;
    });

const Var kRecoShower_EndY(
    [](const caf::SRSliceProxy* slc) -> double {
      double end = -9999.0;
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if (largestShwIdx != -1) {
        double start = slc->reco.shw[largestShwIdx].start.y;
        double dir = slc->reco.shw[largestShwIdx].dir.y;
        double len = slc->reco.shw[largestShwIdx].len;
        end = start + (dir * len);
      }
      return end;
    });

const Var kRecoShower_EndZ(
    [](const caf::SRSliceProxy* slc) -> double {
      double end = -9999.0;
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if (largestShwIdx != -1) {
        double start = slc->reco.shw[largestShwIdx].start.z;
        double dir = slc->reco.shw[largestShwIdx].dir.z;
        double len = slc->reco.shw[largestShwIdx].len;
        end = start + (dir * len);
      }
      return end;
    });

const Var kRecoShower_densityGradient(
    [](const caf::SRSliceProxy* slc) -> double {
      double densityGrad = -5.0;
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if (largestShwIdx != -1) {
        densityGrad = slc->reco.shw[largestShwIdx].selVars.densityGradient;
      }
      return densityGrad;
    });

const Var kRecoShower_densityGradientPower(
    [](const caf::SRSliceProxy* slc) -> double {
      double densityGradPower = -5.0;
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if (largestShwIdx != -1) {
        densityGradPower = slc->reco.shw[largestShwIdx].open_angle;
        densityGradPower = slc->reco.shw[largestShwIdx].selVars.densityGradientPower;
      }
      return densityGradPower;
    });

const Var kRecoShower_trackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      double trackLength = -5.0;
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if (largestShwIdx != -1) {
        trackLength = slc->reco.shw[largestShwIdx].open_angle;
        trackLength = slc->reco.shw[largestShwIdx].selVars.trackLength;
      }
      return trackLength;
    });

const Var kRecoShower_trackWidth(
    [](const caf::SRSliceProxy* slc) -> double {
      double trackWidth = -5.0;
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if (largestShwIdx != -1) {
        trackWidth = slc->reco.shw[largestShwIdx].selVars.trackWidth;
      }
      return trackWidth;
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
      int pdg = -5;
      const int longestTrackIdx(kLongestTrackIdx(slc));
      if (longestTrackIdx != -1) {
        pdg = slc->reco.trk[longestTrackIdx].truth.p.pdg;
      }
      return pdg;
    });

const Var kLongestTrackLength(
    [](const caf::SRSliceProxy* slc) -> double {
      double length = -5.f;
      const int longestTrackIdx(kLongestTrackIdx(slc));
      if (longestTrackIdx != -1) {
        length = slc->reco.trk[longestTrackIdx].len;
      }
      return length;
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
      if (longestTrackIdx != -1 && (kLongestTrackChi2Muon(slc) < 30.f && kLongestTrackChi2Proton(slc) > 60.f)) {
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
