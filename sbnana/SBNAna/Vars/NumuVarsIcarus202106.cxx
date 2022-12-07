#include <cmath>
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana
{
  const Var kTruthIndex = SIMPLEVAR(truth.index);

  const Var kPrimaryEnergy = SIMPLEVAR(truth.E);

  const Var kMuMaxTrack([](const caf::SRSliceProxy* slc) -> float {
      float len(-5.f);
      for (auto const& pfp: slc->reco.pfp)
      {
        const auto& trk = pfp.trk;
	if ( pfp.trackScore > 0.5 && (trk.truth.p.pdg == 13 || trk.truth.p.pdg == -13) && trk.truth.p.length > len ) len = trk.truth.p.length;
      }
      return len;
    });

  const Var kPTrackInd([](const caf::SRSliceProxy* slc) -> int {
      // The (dis)qualification of a slice is based upon the track level information.
      float Longest(0);
      int PTrackInd(-1);
      for (std::size_t i(0); i < slc->reco.pfp.size(); ++i)
      {
	auto const& trk = slc->reco.pfp.at(i).trk;
        if(trk.bestplane == -1 || slc->reco.pfp.at(i).trackScore < 0.5) continue;

	// First we calculate the distance of each track to the slice vertex.
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

	// We require that the distance of the track from the slice is less than
	// 10 cm and that the parent of the track has been marked as the primary.
	const bool AtSlice = ( Atslc < 10.0 && slc->reco.pfp.at(i).parent_is_primary);

        const float Chi2Proton = trk.chi2pid[trk.bestplane].chi2_proton;
        const float Chi2Muon = trk.chi2pid[trk.bestplane].chi2_muon;

	const bool Contained = ( !isnan(trk.end.x) &&
		      ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
		      !isnan(trk.end.y) &&
		      ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		      !isnan(trk.end.z) &&
		      ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
	const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
	const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
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

      if ( kPTrackInd(slc) >= 0 )
      {
	auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
	const bool Contained = ( !isnan(trk.end.x) &&
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
	auto const& trk = slc->reco.pfp.at(kPTrackInd(slc)).trk;
	p = std::hypot(trk.truth.p.genp.x, trk.truth.p.genp.y, trk.truth.p.genp.z);
      }
      return p;
    });


}
