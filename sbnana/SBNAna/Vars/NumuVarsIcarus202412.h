// Stitching function (Reference docdb #38140-v3)

#pragma once

#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"

namespace ana {

    struct PFP {
        // first track
        caf::SRVector3D vertex;      ///< Vertex of slice
        caf::SRVector3D start;       ///< Start point of track
        caf::SRVector3D end;         ///< End point of track
        float len;                   ///< Track length
        float p_muon;                ///< Momentum estimate from trk range (muon hypothesis)
        float slcID;                 ///< Slice ID of the track
        float id;                    ///< PFP ID of the track
        bool muon_1muNp = false;     ///< Muon passing the 1muNp selection (data only)
        bool nuMuCC = false;         ///< Muon from a numuCC interaction (MC only)
        int G4ID;                    ///< G4ID of the track (MC only)
        int pdg;                     ///< PDG of the track (MC only)
        float energy_comp;           ///< Energy completeness of the track (MC only)
        float hit_comp;              ///< Hit completeness of the track (MC only)

        // second track
        caf::SRVector3D start2;      ///< Start point of track
        caf::SRVector3D end2;        ///< End point of track
        float len2;                  ///< Track length
        float p_muon2;               ///< Momentum estimate from trk range (muon hypothesis)
        float slcID2;                ///< Slice ID of the track
        float id2;                   ///< PFP ID of the track
        int G4ID2;                   ///< G4ID of the track (MC only)
        int pdg2;                    ///< PDG of the track (MC only)
        float energy_comp2;          ///< Energy completeness of the track (MC only)
        float hit_comp2;             ///< Hit completeness of the track (MC only)

        // stitching
        size_t nStitch = 0;          ///< Number of stitching in the same event
        float Len;                   ///< Length of the track after the stitching
        float P_muon;                ///< Momentum estimate from the track range (muon hypothesis) after the stitching
    };

    bool kIcarus202401BaryFMCut(const caf::SRSliceProxy &slc) {
        return !std::isnan(slc.barycenterFM.deltaZ_Trigger) && 
            slc.barycenterFM.deltaZ_Trigger >= 0 && 
            slc.barycenterFM.deltaZ_Trigger < 100;
    }

    bool kIcarus202412RecoFiducial(const caf::SRSliceProxy &slc) {
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
    bool kIcarus202401MuonTrack(const caf::SRPFPProxy &pfp, float trkScore, float trkLmin) {
        // The (dis)qualification of a slice is based upon the track level information.
        bool muTrack = false;
        if(pfp.trackScore < trkScore) return muTrack;
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
        const bool MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len >= trkLmin );
        if(MaybeMuonContained) muTrack = true;
        return muTrack;
    }

    extern const SpillMultiVar kIcarus202412Stitch;

}