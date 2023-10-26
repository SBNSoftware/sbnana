

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"

//#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Vars/NumuVars.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/NueCuts.h"

#include "TMath.h"
//#include "TVector3.h"
#include "TRotation.h"

using namespace ana;

// A Var is a little snippet of code that takes a "proxy" representing the
// event record and returns a single number to plot

const caf::SRPFPProxy* LargestRecoShower(const caf::SRSliceProxy* slc){  
  double largestE = 0;
  const caf::SRPFPProxy* largestShower;
  for (const auto& pfp : slc->reco.pfp) {
    if (pfp.trackScore > 0.5) { continue; }
    if (pfp.shw.bestplane_energy > largestE) {
      largestE = pfp.shw.bestplane_energy;
      largestShower= &pfp;
    }
  }
  if (largestE == 0) { return 0; }
  return largestShower;
}

const MultiVar kMatches([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> matches;
  const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
  if (!largestShower) { return {}; }
  for (const auto& pfp : slc->reco.pfp) {
    if (pfp.trackScore > 0.5) { continue; }
    if (pfp.id == largestShower->id) { continue; }
    if (pfp.shw.bestplane_energy < 0) { continue; } 
    bool isMatch = (pfp.shw.truth.bestmatch.G4ID == largestShower->shw.truth.bestmatch.G4ID);
    if (isMatch) { matches.push_back(1); }
    else { matches.push_back(0); }
  }
  return matches;
});

const MultiVar kEnergies([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> energies;
  const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
  if (!largestShower) { return {}; }
  for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
    if (pfp.id == largestShower->id) { continue; }
    if (pfp.shw.bestplane_energy < 0) { continue; } 
    energies.push_back(pfp.shw.bestplane_energy);
  }
  return energies;
});

const MultiVar kSlices([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> slices;
  const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
  if (!largestShower) { return {}; }
  for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
    if (pfp.id == largestShower->id) { continue; }
    if (pfp.shw.bestplane_energy < 0) { continue; } 
    if (!slc->is_clear_cosmic) { slices.push_back(1); }
    else { slices.push_back(0); }
  }
  return slices;
});

const MultiVar kAngles([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> angles;
  const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
  if (!largestShower) { return {}; }
  for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
      if (pfp.id == largestShower->id) { continue; }
      if (pfp.shw.bestplane_energy < 0) { continue; } 
    TVector3 recoDir(pfp.shw.end.x-pfp.shw.start.x,
 		 pfp.shw.end.y-pfp.shw.start.y,
		 pfp.shw.end.z-pfp.shw.start.z);
    recoDir = recoDir.Unit();
    TVector3 largestDir(largestShower->shw.end.x-largestShower->shw.start.x,
		 largestShower->shw.end.y-largestShower->shw.start.y,
		 largestShower->shw.end.z-largestShower->shw.start.z);
    largestDir = largestDir.Unit();
    double dot = largestDir.Dot(recoDir);
    double theta = TMath::ACos(dot);
    angles.push_back(theta);
  }
  return angles;
});

const MultiVar kEndtoend([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> dist;
    const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
    if (!largestShower) { return {}; }
    for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
      if (pfp.id == largestShower->id) { continue; }
      if (pfp.shw.bestplane_energy < 0) { continue; } 
      TVector3 recoTraj(largestShower->shw.end.x-pfp.shw.end.x,
  	  	        largestShower->shw.end.y-pfp.shw.end.y,
		        largestShower->shw.end.z-pfp.shw.end.z);
      double distance = std::abs(recoTraj.Mag());
      dist.push_back(distance);
    }
  return dist;
});

const MultiVar kEndtostart([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> dist;
    const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
    if (!largestShower) { return {}; }
    for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
      if (pfp.id == largestShower->id) { continue; }
      if (pfp.shw.bestplane_energy < 0) { continue; } 
      TVector3 recoTraj(largestShower->shw.end.x-pfp.shw.start.x,
  	  	        largestShower->shw.end.y-pfp.shw.start.y,
		        largestShower->shw.end.z-pfp.shw.start.z);
      double distance = std::abs(recoTraj.Mag());
      dist.push_back(distance);
    }
  return dist;
});

const MultiVar kStarttostart([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> dist;
    const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
    if (!largestShower) { return {}; }
    for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
      if (pfp.id == largestShower->id) { continue; }
      if (pfp.shw.bestplane_energy < 0) { continue; } 
      TVector3 recoTraj(largestShower->shw.start.x-pfp.shw.start.x,
  	  	        largestShower->shw.start.y-pfp.shw.start.y,
		        largestShower->shw.start.z-pfp.shw.start.z);
      double distance = std::abs(recoTraj.Mag());
      dist.push_back(distance);
    }
  return dist;
});

const MultiVar kStarttoend([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> dist;
    const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
    if (!largestShower) { return {}; }
    for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
      if (pfp.id == largestShower->id) { continue; }
      if (pfp.shw.bestplane_energy < 0) { continue; } 
      TVector3 recoTraj(largestShower->shw.start.x-pfp.shw.end.x,
  	  	        largestShower->shw.start.y-pfp.shw.end.y,
		        largestShower->shw.start.z-pfp.shw.end.z);
      double distance = std::abs(recoTraj.Mag());
      dist.push_back(distance);
    }
  return dist;
});

const MultiVar kIntersections([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> intersections;
    const caf::SRPFPProxy* largestShower = LargestRecoShower(slc);
    if (!largestShower) { return {}; }
    for (const auto& pfp : slc->reco.pfp) {
      if (pfp.trackScore > 0.5) { continue; }
      if (pfp.id == largestShower->id) { continue; }
      if (pfp.shw.bestplane_energy <= 0) { continue; } 

      TVector3 r(largestShower->shw.start.x-pfp.shw.start.x,
  	  	        largestShower->shw.start.y-pfp.shw.start.y,
		        largestShower->shw.start.z-pfp.shw.start.z);
      TVector3 r_l(pfp.shw.start.x-pfp.shw.end.x,
  	  	        pfp.shw.start.y-pfp.shw.end.y,
		        pfp.shw.start.z-pfp.shw.end.z);
      TVector3 r_L(largestShower->shw.start.x-largestShower->shw.end.x,
  	  	        largestShower->shw.start.y-largestShower->shw.end.y,
		        largestShower->shw.start.z-largestShower->shw.end.z);

      TVector3 orth(r_l.Cross(r_L));
      orth = orth.Unit();

      double theta_L = largestShower->shw.open_angle;
      double theta_l = pfp.shw.open_angle;

      double length_L = largestShower->shw.len/TMath::Cos(theta_L);
      double length_l = pfp.shw.len/TMath::Cos(theta_l);

      TRotation rot_L1;
      TRotation rot_L2;
      TRotation rot_l1;
      TRotation rot_l2;

      rot_L1.Rotate(theta_L,orth);
      rot_L2.Rotate(2*TMath::Pi()-theta_L,orth);
      rot_l1.Rotate(theta_l,orth);
      rot_l2.Rotate(2*TMath::Pi()-theta_l,orth);

      TVector3 r_L1 = (rot_L1*r_L);
      TVector3 r_L2 = (rot_L2*r_L);
      TVector3 r_l1 = (rot_l1*r_l);
      TVector3 r_l2 = (rot_l2*r_l);

      r_L1.SetMag(2*length_L); 
      r_L2.SetMag(2*length_L); 
      r_l1.SetMag(2*length_l); 
      r_l2.SetMag(2*length_l);

      TVector3 P1(pfp.shw.start.x, pfp.shw.start.y, pfp.shw.start.z);
      TVector3 Q1(largestShower->shw.start.x, largestShower->shw.start.y, largestShower->shw.start.z);  

      std::vector<TVector3> smallShws({r_l1,r_l2});
      std::vector<TVector3> largeShws({r_L1,r_L2});

      std::vector<int> ints;
      int counter = 0;
      for (const auto& small_cone : smallShws) {
        TVector3 P2 = P1 + small_cone;
        for (const auto& large_cone : largeShws) {
          TVector3 Q2 = Q1 + large_cone;
 
          double a = (P2 - P1).Dot(Q1-P1)/((P2-P1).Mag() * (P2-P1).Mag());
          double b = (P2 - P1).Dot(Q2-Q1)/((P2-P1).Mag() * (P2-P1).Mag());
          TVector3 c = b * (P2-P1)-(Q2-Q1);

          double t1 = c.Dot(Q1-(1-a)*P1 - a*P2)/(c.Mag() * c.Mag());
          double t0 = a + t1*b;
          if (t0 > 0 && t0 < 1 && t1 > 0 && t1 < 1) { 
            ints.push_back(1);
            counter++;
          }
          else { ints.push_back(0); }
        }
      }

      intersections.push_back(counter);
    }
  return intersections;
});


