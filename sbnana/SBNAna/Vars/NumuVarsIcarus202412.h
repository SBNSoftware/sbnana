// Stitching function (Reference docdb #38140-v3)

#pragma once

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"

namespace ana {

    struct PFP {
        caf::SRVector3D vertex;  ///< Vertex of slice
        caf::SRVector3D start;   ///< Start point of track
        caf::SRVector3D end;     ///< End point of track

        float len;
        float slcID;
        float id;

        float G4ID;
        float pdg;
        float energy_comp;
        float hit_comp;
    };

    bool kIcarus202401RecoFiducial(const caf::SRSliceProxy &slc) {
        return ( !isnan(slc.vertex.x) &&
	        ( ( slc.vertex.x < -61.94 - 25 && slc.vertex.x > -358.49 + 25 ) ||
		      ( slc.vertex.x > 61.94 + 25 && slc.vertex.x < 358.49 - 25 ) ) &&
	        !isnan(slc.vertex.y) &&
	        ( slc.vertex.y > -181.86 + 25 && slc.vertex.y < 134.96 - 25 ) &&
	        !isnan(slc.vertex.z) &&
	        ( slc.vertex.z > -894.95 + 30 && slc.vertex.z < 894.95 - 50 ) &&
            !(slc.vertex.x > 210 && slc.vertex.y > 60 && slc.vertex.z > 290 && slc.vertex.z < 390)
        );
    }

    // adapted from kIcarus202401MuonIdx in NumuVarsIcarus202401.cxx
    bool kIcarus202401MuonTrack(const caf::SRPFPProxy &pfp) {
        // The (dis)qualification of a slice is based upon the track level information.
        bool muTrack = false;
        if(pfp.trackScore < 0.4) return muTrack;
        if(std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x)) return muTrack;
        auto const& trk = pfp.trk;

        //int plane = trk.calo[1].nhit > trk.calo[2].nhit ? 1 : 2;
        int plane = 2; // Hard code collection plane for now since induction 2 has peak at higher chi2

        //float Chi2Proton = trk.chi2pid[plane].chi2_proton;
        //float Chi2Muon = trk.chi2pid[plane].chi2_muon;
        auto chi2 = chi2pid::chi2_calculator.calculate_chi2(trk.calo[plane]);
        float Chi2Proton = chi2.chi2_proton;
        float Chi2Muon = chi2.chi2_muon;

        const bool Contained = ( !isnan(trk.end.x) &&
            (( trk.end.x < -61.94 - 5 && trk.end.x > -358.49 + 5 ) ||
             ( trk.end.x > 61.94 + 5 && trk.end.x < +358.49 - 5 )) &&
            !isnan(trk.end.y) &&
            ( trk.end.y > -181.86 + 5 && trk.end.y < 134.96 - 5 ) &&
            !isnan(trk.end.z) &&
            ( trk.end.z > -894.95 + 5 && trk.end.z < 894.95 - 5 ) );
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len >= 20 );
        if(MaybeMuonContained) muTrack = true;
        return muTrack;
    }

    extern const SpillMultiVar kIcarus202412Stitch;

}