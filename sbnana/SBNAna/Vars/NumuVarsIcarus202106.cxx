#include <cmath>
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  const Var kTruthIndex = SIMPLEVAR(truth.index);
  
  const Var kPrimaryEnergy = SIMPLEVAR(truth.E);
  
  const Var kMuMaxTrack([](const caf::SRSliceProxy* slc) -> float {
      float len(-5.f);
      for (auto const& trk : slc->reco.trk)
      {
	if ( (trk.truth.p.pdg == 13 || trk.truth.p.pdg == -13) && trk.truth.p.length > len ) len = trk.truth.p.length;
      }
      return len;
    });

  const Var kPTrackInd([](const caf::SRSliceProxy* slc) -> int {
      // The (dis)qualification of a slice is based upon the track level information.
      float Atslc, Chi2Proton, Chi2Muon, Longest(0);
      bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained;
      int PTrackInd(-1);
      for (std::size_t i(0); i < slc->reco.trk.size(); ++i)
      {
	auto const& trk = slc->reco.trk.at(i);
	// First we calculate the distance of each track to the slice vertex.
	Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
		      pow( slc->vertex.y - trk.start.y, 2.0 ) +
		      pow( slc->vertex.z - trk.start.z, 2.0 ) );
	// We require that the distance of the track from the slice is less than
	// 10 cm and that the parent of the track has been marked as the primary.
	AtSlice = ( Atslc < 10.0);// && trk.pfp.parent_is_primary);
	
	if (trk.bestplane == 0)
	{
	  Chi2Proton = trk.chi2pid0.chi2_proton;
	  Chi2Muon = trk.chi2pid0.chi2_muon;
	}
	else if (trk.bestplane == 1)
	{
	  Chi2Proton = trk.chi2pid1.chi2_proton;
	  Chi2Muon = trk.chi2pid1.chi2_muon;
	}
	else
	{
	  Chi2Proton = trk.chi2pid2.chi2_proton;
	  Chi2Muon = trk.chi2pid2.chi2_muon;
	}

	Contained = ( !isnan(trk.end.x) &&
		      ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
		      !isnan(trk.end.y) &&
		      ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		      !isnan(trk.end.z) &&
		      ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
	MaybeMuonExiting = ( !Contained && trk.len > 100);
	MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
	if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
	{
	  Longest = trk.len;
	  PTrackInd = i;
	}
      }
      return PTrackInd;
    });

  const Var kRecoMuonP([](const caf::SRSliceProxy* slc) -> float {
      float p(-5.f);
      bool Contained(false);
      
      if ( kPTrackInd(slc) >= 0 )
      {
	auto const& trk = slc->reco.trk.at(kPTrackInd(slc));
	Contained = ( !isnan(trk.end.x) &&
		      ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
		      !isnan(trk.end.y) &&
		      ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		      !isnan(trk.end.z) &&
		      ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
	if(Contained) p = trk.rangeP.p_muon;
	else p = trk.mcsP.fwdP_muon;
      }
      return p;
    });

  const Var kTrueMuonP([](const caf::SRSliceProxy* slc) -> float {
      float p(-5.f);
      
      if ( kPTrackInd(slc) >= 0 )
      {
	auto const& trk = slc->reco.trk.at(kPTrackInd(slc));
	p = sqrt( pow(trk.truth.p.genp.x, 2) + 
		  pow(trk.truth.p.genp.y, 2) + 
		  pow(trk.truth.p.genp.z, 2) );
      }
      return p;
    });


}
