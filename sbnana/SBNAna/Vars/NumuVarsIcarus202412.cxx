// Stitching function (Reference docdb #38140-v3)

#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202412.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"

// ROOT
#include "TFile.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TProfile.h"

// C++
#include <vector>
#include <algorithm>
#include <iostream>

namespace ana {

    namespace slc {
        // parameters
        bool print = false;  // print on terminal
        bool EHcomp = false;  // enable quality cut on hit_comp and E_comp
        bool DvCut = false;  // enable cut on the distance between mc vtx and reco vtx
        float mcLmin = 50.;  // minimal length of the track on mc
        float Vdist = 15.;   // distance between mc vtx and reco vtx
        float Tss = 0.7;     // overlap of the segments on the z-axis
        float thetaM = 35.;  // angle between the directions of the segments
        float Dss = 0.7;     // 3D distance between the segments
    } // namespace slc

    struct PFP {
        caf::SRVector3D start;       ///< Start point of track
        caf::SRVector3D end;         ///< End point of track

        float len;
        float slcID;
        float id;

        float G4ID;
        float pdg;
        float energy_comp;
        float hit_comp;
    };

    static bool kIcarus202401BaryFMCut(const caf::SRSliceProxy &slc) {
        return !std::isnan(slc.barycenterFM.deltaZ_Trigger) && 
            slc.barycenterFM.deltaZ_Trigger >= 0 && 
            slc.barycenterFM.deltaZ_Trigger < 100;
    }

    static bool kIcarus202401RecoFiducial(const caf::SRSliceProxy &slc) {
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

    // sort start and end points of the track according to the distance from the vertex
    void Icarus202412SortPFP(PFP &a, const caf::SRVector3DProxy &vtx) {
        float dist[2];
        dist[0] = std::hypot(a.start.x-vtx.x, a.start.y-vtx.y, a.start.z-vtx.z);
        dist[1] = std::hypot(a.end.x-vtx.x, a.end.y-vtx.y, a.end.z-vtx.z);
        if(dist[1] < dist[0]) std::swap(a.start, a.end);
    }

    // stitching function within the slice
    std::vector<double> Icarus202412Stitch(const caf::SRSliceProxy &slc, std::vector<double> &hist, std::vector<int> &sthG4ID, bool mc) {
        std::vector<double> stitch;
        bool found = false;
        bool wrong = false;
        int muG4ID = 0;
        size_t nStitch = 0;
        // check if the slice is associated with an OptFlash
        if(!kIcarus202401BaryFMCut(slc)) return stitch;
        // check if the slice is in the fiducial volume
        if(!kIcarus202401RecoFiducial(slc)) return stitch;
        auto const& vtx = slc.vertex;
        // MC
        if(mc) {
            // check on the distance between mc vtx and reco vtx
            if(slc::DvCut && !std::isnan(slc.truth.position.x)) {
                auto const& mcvtx = slc.truth.position;
                float dist = std::hypot(mcvtx.x-vtx.x, mcvtx.y-vtx.y, mcvtx.z-vtx.z);
                if(dist > slc::Vdist) return stitch;
            }
            // select mu primary daugthers of numuCC interactions, contained in a single cryostat
            // and with a minimal length of the track
            for(int i=0; i<slc.truth.nprim; i++) {
                if(slc.truth.prim.at(i).pdg != 13 && slc.truth.prim.at(i).pdg != -13) continue;
                if(!slc.truth.prim.at(i).contained) continue;
                if(slc.truth.prim.at(i).length < slc::mcLmin) continue;
                muG4ID = slc.truth.prim.at(i).G4ID;
            }
        }
        size_t counter = 0;
        size_t nsegments = 0;
        // loop on the pfp
        std::vector<PFP> pfp;
        for(size_t i=0; i<slc.reco.npfp; i++) {
            // check if the pfp is a muon
            if(!kIcarus202401MuonTrack(slc.reco.pfp.at(i))) continue;
            PFP P;
            P.start.x = slc.reco.pfp.at(i).trk.start.x;
            P.start.y = slc.reco.pfp.at(i).trk.start.y;
            P.start.z = slc.reco.pfp.at(i).trk.start.z;
            P.end.x = slc.reco.pfp.at(i).trk.end.x;
            P.end.y = slc.reco.pfp.at(i).trk.end.y;
            P.end.z = slc.reco.pfp.at(i).trk.end.z;
            P.len = slc.reco.pfp.at(i).trk.len;
            P.slcID = slc.reco.pfp.at(i).slcID;
            P.id = slc.reco.pfp.at(i).id;
            // MC
            if(mc) {
                P.G4ID = slc.reco.pfp.at(i).trk.truth.bestmatch.G4ID;
                P.pdg = slc.reco.pfp.at(i).trk.truth.p.pdg;
                P.energy_comp = slc.reco.pfp.at(i).trk.truth.bestmatch.energy_completeness;
                P.hit_comp = slc.reco.pfp.at(i).trk.truth.bestmatch.hit_completeness;
                if(muG4ID == P.G4ID) {
                    if(P.energy_comp >= 0.1 && P.hit_comp >= 0.1) counter++;
                    nsegments++;
                }
            }
            pfp.push_back(P);
        }
        // quality cut on hit_comp and E_comp
        if(slc::EHcomp && (counter < 2 || nsegments != 2)) return stitch;
        // order the pfp according to the length
        std::sort(pfp.begin(), pfp.end(), [](const PFP a, const PFP b) {return a.len > b.len;});
        // stitching
        for(size_t i=0; i<pfp.size(); i++) {
            for(size_t j=i+1; j<pfp.size(); j++) {
                std::vector<PFP> trk;
                // sort vertex and segments according to the baycenters
                std::vector<float> ordz;
                ordz.push_back(vtx.z);                                       // vertex
                ordz.push_back( 0.5*(pfp.at(i).start.z+pfp.at(i).end.z) );   // (start 1 + end 1)/2
                ordz.push_back( 0.5*(pfp.at(j).start.z+pfp.at(j).end.z) );   // (start 2 + end 2)/2
                std::sort(ordz.begin(), ordz.end());
                // check vertex in the middle
                if(ordz[1] == vtx.z) {
                    if(mc) {
                        if(pfp.at(i).G4ID==muG4ID && pfp.at(j).G4ID==muG4ID && found == false) {
                            hist.push_back(1);
                            hist.push_back(4);
                            found = true;
                        }
                        else if(pfp.at(i).G4ID==muG4ID || pfp.at(j).G4ID==muG4ID) hist.push_back(9);
                    }
                    continue;
                }
                // check which segment is closest to the vertex
                if(ordz[1]==0.5*(pfp.at(i).start.z+pfp.at(i).end.z)) {
                    trk.push_back(pfp.at(i));
                    trk.push_back(pfp.at(j));
                }
                else {
                    trk.push_back(pfp.at(j));
                    trk.push_back(pfp.at(i));
                }
                // sort start and end points of each segment according to distance from the vertex
                for(size_t k=0; k<2; k++) Icarus202412SortPFP(trk.at(k), vtx);
                // check segment-segment overlap on z-axis
                std::vector<float> zlen;
                zlen.push_back(std::abs(trk.at(0).end.z-trk.at(0).start.z));
                zlen.push_back(std::abs(trk.at(1).end.z-trk.at(1).start.z));
                std::sort(zlen.begin(),zlen.end());
                bool Zoverlap = (ordz[0]==vtx.z && trk.at(1).start.z < trk.at(0).end.z && std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= slc::Tss) ||
                                (ordz[2]==vtx.z && trk.at(1).start.z > trk.at(0).end.z && std::abs(trk.at(1).start.z - trk.at(0).end.z)/zlen[0] >= slc::Tss);
                if(Zoverlap) {
                    if(mc) {
                        if(trk.at(0).G4ID==muG4ID && trk.at(1).G4ID==muG4ID && found == false) {
                            hist.push_back(1);
                            hist.push_back(5);
                            found = true;
                        }
                        else if(trk.at(0).G4ID==muG4ID || trk.at(1).G4ID==muG4ID) hist.push_back(10);
                    }
                    continue;
                }
                // check angle between segments
                float a[3],b[3];
                a[0] = trk.at(0).end.x-trk.at(0).start.x;
                a[1] = trk.at(0).end.y-trk.at(0).start.y;
                a[2] = trk.at(0).end.z-trk.at(0).start.z;
                b[0] = trk.at(1).end.x-trk.at(1).start.x;
                b[1] = trk.at(1).end.y-trk.at(1).start.y;
                b[2] = trk.at(1).end.z-trk.at(1).start.z;
                float Ra = std::hypot(a[0],a[1],a[2]);
                float Rb = std::hypot(b[0],b[1],b[2]);
                float c = (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/(Ra*Rb); // scalar product
                float phi = acos(c); // [rad]
                if(phi*180/3.141592 >= slc::thetaM) {
                    if(mc) {
                        if(trk.at(0).G4ID==muG4ID && trk.at(1).G4ID==muG4ID && found == false) {
                            hist.push_back(1);
                            hist.push_back(6);
                            found = true;
                        }
                        else if(trk.at(0).G4ID==muG4ID || trk.at(1).G4ID==muG4ID) hist.push_back(11);
                    }
                    continue;
                }
                // check distance between segments
                std::vector<float> length;
                length.push_back(Ra);
                length.push_back(Rb);
                std::sort(length.begin(), length.end());
                float dist = std::hypot(trk.at(1).start.x-trk.at(0).end.x, trk.at(1).start.y-trk.at(0).end.y, trk.at(1).start.z-trk.at(0).end.z);
                if(dist/length[0] >= slc::Dss) {
                    if(mc) {
                        if(trk.at(0).G4ID==muG4ID && trk.at(1).G4ID==muG4ID && found == false) {
                            hist.push_back(1);
                            hist.push_back(7);
                            found = true;
                        }
                        else if(trk.at(0).G4ID==muG4ID || trk.at(1).G4ID==muG4ID) hist.push_back(12);
                    }
                    continue;
                }
                // fill last bins
                if(mc) {
                    if(trk.at(0).G4ID==muG4ID && trk.at(1).G4ID==muG4ID && found==false) {
                        hist.push_back(0);
                        found = true;
                    }
                    else if(!(trk.at(0).G4ID==muG4ID && trk.at(1).G4ID==muG4ID) && (trk.at(0).G4ID==muG4ID || trk.at(1).G4ID==muG4ID) && wrong==false) {
                        hist.push_back(2);
                        wrong = true;
                    }
                }
                nStitch++;
                stitch.push_back(nStitch);                           // number of stitching
                stitch.push_back(trk.at(0).slcID);                   // slcID
                stitch.push_back(trk.at(0).id);                      // pfpID 1
                stitch.push_back(trk.at(1).id);                      // pfpID 2

                stitch.push_back(vtx.x);                             // vertex
                stitch.push_back(trk.at(0).start.x);                 // start 1
                stitch.push_back(trk.at(0).end.x);                   // end 1
                stitch.push_back(trk.at(1).start.x);                 // start 2
                stitch.push_back(trk.at(1).end.x);                   // end 2

                stitch.push_back(vtx.y);                             // vertex
                stitch.push_back(trk.at(0).start.y);                 // start 1
                stitch.push_back(trk.at(0).end.y);                   // end 1
                stitch.push_back(trk.at(1).start.y);                 // start 2
                stitch.push_back(trk.at(1).end.y);                   // end 2

                stitch.push_back(vtx.z);                             // vertex
                stitch.push_back(trk.at(0).start.z);                 // start 1
                stitch.push_back(trk.at(0).end.z);                   // end 1
                stitch.push_back(trk.at(1).start.z);                 // start 2
                stitch.push_back(trk.at(1).end.z);                   // end 2

                stitch.push_back(trk.at(0).len);                     // length before stitching
                stitch.push_back(trk.at(0).len + trk.at(1).len);     // length after stitching
                if(mc) sthG4ID.push_back(muG4ID);                    // mu G4ID
            }
        }
        if(found == true) hist.push_back(14);
        return stitch;
    }

    // print on terminal the output of stitch
    void Icarus202412PrintStitch(std::vector<double> &stitch) {
        for(size_t i=0; i<stitch.size(); i=i+21) {
            std::cout << " " << std::endl;
            std::cout << "===================" << std::endl;
            std::cout << "Stitching No: " << stitch.at(0+i) << std::endl;
            std::cout << " " << std::endl;
            std::cout << "slcID = " << stitch.at(1+i) << std::endl;
            std::cout << "pfpID 1 = " << stitch.at(2+i) << std::endl;
            std::cout << "pfpID 2 = " << stitch.at(3+i) << std::endl;
            std::cout << " " << std::endl;
            std::cout << "geom = {v,s1,e1,s2,e2}" << std::endl;
            for(size_t j=0; j<5; j++) {
                std::cout << "x[" << j << "] = " << stitch.at(4+j+i) <<
                    "\t y[" << j << "] = " << stitch.at(9+j+i) <<
                    "\t z[" << j << "] = " << stitch.at(14+j+i) << std::endl;
            }
            std::cout << " " << std::endl;
            std::cout << "L1 = " << stitch.at(19+i) << std::endl;
            std::cout << "L1+L2 = " << stitch.at(20+i) << std::endl;
            std::cout << " " << std::endl;
        }
    }

    // stitching function within the slice for sbnana
    const SpillMultiVar kIcarus202412Stitch([](const caf::SRSpillProxy* sr) -> std::vector<double> {
        std::vector<double> stitch;
        std::vector<double> hist;
        std::vector<int> sthG4ID;
        bool mc = false;
        if(sr->hdr.ismc) mc = true;
        for(int i=0; i<sr->nslc; i++) {
            std::vector<double> temp;
            temp = Icarus202412Stitch(sr->slc.at(i), hist, sthG4ID, mc);
            if(!temp.empty() && slc::print) {
                std::cout << " " << std::endl;
                std::cout << " " << std::endl;
                std::cout << "==============" << std::endl;
                std::cout << "Run: " << sr->hdr.run << std::endl;
                std::cout << "Subrun: " << sr->hdr.subrun << std::endl;
                std::cout << "Event: " << sr->hdr.evt << std::endl;
                std::cout << "==============" << std::endl;
                Icarus202412PrintStitch(temp);
            }
            for(size_t j=0; j<temp.size(); j++) stitch.push_back(temp.at(j));
        }
        return hist;
    });

} // namespace ana