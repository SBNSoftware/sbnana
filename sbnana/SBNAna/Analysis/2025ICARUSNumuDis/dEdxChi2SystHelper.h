#pragma once

// Helper for the dE/dx wireMod like systematic.
// Artificially shifts dE/dx hit values and recomputes chi2 PID scores.

#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "cetlib/search_path.h"

#include "TMath.h"
#include "TFile.h"
#include "TProfile.h"

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>
#include <algorithm>

namespace ana {

namespace dedxsyst {


// Modified Box Model parameters from icaruscode/fcl/reco/Definitions/stage1_icarus_defs.fcl
// (G. Putnam, Nov 2022 -- docdb 28639)
// ModBoxA = 0.904
// ModBoxBTF1: [0]/sqrt(sin^2(phi) + cos^2(phi)/[1]^2)  with [0]=0.204, [1]=1.25
static constexpr double BoxA  = 0.904;
static constexpr double BoxB0 = 0.204;
static constexpr double BoxBk = 1.25;

inline double getBeta(double phi)
{
    double phi_rad = phi * TMath::Pi() / 180.;
    return BoxB0 / std::sqrt(std::pow(std::sin(phi_rad), 2) +
                              std::pow(std::cos(phi_rad), 2) / (BoxBk * BoxBk));
}

inline double CorrecteddEdx(double dEdx, double dirx, double f)
{
    double phi  = std::acos(std::abs(dirx)) * 180. / TMath::Pi();
    double beta = getBeta(phi);
    return (std::pow(beta * dEdx + BoxA, f) - BoxA) / beta;
}

struct Chi2PIDCorrected {
    Chi2PIDCorrected()
        : file_Chi2Template(nullptr),
          dedx_range_pro(nullptr), dedx_range_ka(nullptr),
          dedx_range_pi(nullptr),  dedx_range_mu(nullptr) {}

    ~Chi2PIDCorrected() {
        if (file_Chi2Template) file_Chi2Template->Close();
    }

    void check_graphs() {
        if (dedx_range_pro) return;
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fname = "dEdxrestemplates.root";
        std::string fpath;
        if (!sp.find_file(fname, fpath)) {
            std::cerr << "dEdxChi2SystHelper: cannot find '" << fname << "'\n";
            std::abort();
        }
        file_Chi2Template = TFile::Open(fpath.c_str());
        dedx_range_pro = file_Chi2Template->Get<TProfile>("dedx_range_pro");
        dedx_range_ka  = file_Chi2Template->Get<TProfile>("dedx_range_ka");
        dedx_range_pi  = file_Chi2Template->Get<TProfile>("dedx_range_pi");
        dedx_range_mu  = file_Chi2Template->Get<TProfile>("dedx_range_mu");
    }

    chi2pid::Chi2PID::Chi2PIDResult
    calculate_chi2_corrected(const caf::Proxy<caf::SRTrackCalo>& calo,
                             double dirx, double f)
    {
        check_graphs();

        chi2pid::Chi2PID::Chi2PIDResult output{0, 0, 0, 0, 0, 0};
        if (calo.points.size() == 0) return output;

        int npt = 0;
        double chi2pro = 0, chi2ka = 0, chi2pi = 0, chi2mu = 0;
        double PIDA = 0;
        std::vector<double> vpida;
        int used_trkres = 0;

        for (unsigned i = 1; i < calo.points.size() - 1; ++i) {
            const auto& pt = calo.points[i];
            double hit_dedx = CorrecteddEdx(pt.dedx, dirx, f);
            double hit_rr   = pt.rr;

            if (hit_rr > 25) continue;
            if (hit_rr < 25) {
                PIDA += hit_dedx * pow(hit_rr, 0.42);
                vpida.push_back(hit_dedx * pow(hit_rr, 0.42));
                used_trkres++;
            }
            if (hit_dedx > 100 || hit_dedx < 0.5) continue;

            int bin = dedx_range_pro->FindBin(hit_rr);
            if (bin < 1 || bin > dedx_range_pro->GetNbinsX()) continue;

            auto bin_content = [&](TProfile* h) {
                double c = h->GetBinContent(bin);
                if (c < 1e-6)
                    c = (h->GetBinContent(bin-1) + h->GetBinContent(bin+1)) / 2.;
                return c;
            };
            auto bin_error = [&](TProfile* h) {
                double e = h->GetBinError(bin);
                if (e < 1e-6)
                    e = (h->GetBinError(bin-1) + h->GetBinError(bin+1)) / 2.;
                return e;
            };

            double bincpro = bin_content(dedx_range_pro);
            double bincka  = bin_content(dedx_range_ka);
            double bincpi  = bin_content(dedx_range_pi);
            double bincmu  = bin_content(dedx_range_mu);

            double binepro = bin_error(dedx_range_pro);
            double bineka  = bin_error(dedx_range_ka);
            double binepi  = bin_error(dedx_range_pi);
            double binemu  = bin_error(dedx_range_mu);

            double errdedx = (0.04231 + 0.0001783 * hit_dedx * hit_dedx) * hit_dedx;

            auto chi2term = [&](double binc, double bine) {
                return pow((hit_dedx - binc) / std::sqrt(pow(bine,2) + pow(errdedx,2)), 2);
            };
            chi2pro += chi2term(bincpro, binepro);
            chi2ka  += chi2term(bincka,  bineka);
            chi2pi  += chi2term(bincpi,  binepi);
            chi2mu  += chi2term(bincmu,  binemu);
            ++npt;
        }

        if (npt) {
            chi2pro /= npt; chi2ka /= npt;
            chi2pi  /= npt; chi2mu /= npt;
        }
        if (used_trkres > 0)
            PIDA = TMath::Median(vpida.size(), &vpida[0]);

        output.chi2_proton = chi2pro;
        output.chi2_kaon   = chi2ka;
        output.chi2_pion   = chi2pi;
        output.chi2_muon   = chi2mu;
        output.PIDA        = PIDA;
        output.ndof        = npt;
        return output;
    }

    TFile*    file_Chi2Template;
    TProfile *dedx_range_pro, *dedx_range_ka, *dedx_range_pi, *dedx_range_mu;
};

inline Chi2PIDCorrected chi2_corrected_calculator;

// --------------------------------------------------------------------------
// WireMod scale factor lookup
// Files: SBNAna/WireMod/ICARUS/Run{2,4}/XW_Coarse_TPC{0-3}_plane{0-2}.txt
// Columns: theta_low theta_high x_low x_high theta_center x_center
//          ratio_integral err_integral ratio_width err_width
// --------------------------------------------------------------------------

// Navigate from this header's location (.../SBNAna/Analysis/2025ICARUSNumuDis/)
// up to SBNAna/ and then into WireMod/ICARUS/
inline std::string wiremod_base_dir()
{
    std::string f = __FILE__;
    for (int i = 0; i < 3; ++i) f = f.substr(0, f.rfind('/'));
    return f + "/WireMod/ICARUS/";
}

struct WireModBin {
    float theta_low, theta_high;
    float x_low, x_high;
    float ratio_integral, err_integral;
};

struct WireModTable {
    // (run_period, tpc, plane) -> bins
    std::map<std::tuple<int,int,int>, std::vector<WireModBin>> data;
    std::set<int> loaded;

    void load(int run_period)
    {
        if (loaded.count(run_period)) return;
        loaded.insert(run_period);
        std::string base = wiremod_base_dir() + "Run" + std::to_string(run_period) + "/";
        for (int tpc = 0; tpc < 4; ++tpc) {
            for (int plane = 0; plane < 3; ++plane) {
                std::string fname = base + "XW_Coarse_TPC" + std::to_string(tpc)
                                        + "_plane" + std::to_string(plane) + ".txt";
                std::ifstream fin(fname);
                if (!fin) {
                    std::cerr << "WireModTable: cannot open " << fname << "\n";
                    continue;
                }
                auto& bins = data[{run_period, tpc, plane}];
                std::string line;
                while (std::getline(fin, line)) {
                    if (line.empty() || line[0] == '#') continue;
                    std::istringstream ss(line);
                    WireModBin b;
                    float theta_center, x_center, ratio_width, err_width;
                    ss >> b.theta_low >> b.theta_high
                       >> b.x_low    >> b.x_high
                       >> theta_center >> x_center
                       >> b.ratio_integral >> b.err_integral
                       >> ratio_width >> err_width;
                    bins.push_back(b);
                }
            }
        }
    }

    // Returns ratio_integral for the bin matching (thetaXW, x), or 1.0 if out of range.
    double scale_factor(double thetaXW, double x, int tpc, int plane, int run_period)
    {
        load(run_period);
        auto it = data.find({run_period, tpc, plane});
        if (it == data.end()) return 1.0;
        for (const auto& b : it->second) {
            if (thetaXW >= b.theta_low && thetaXW < b.theta_high &&
                x       >= b.x_low     && x       < b.x_high)
                return b.ratio_integral;
        }
        return 1.0;
    }
};

inline WireModTable wiremod_table;

// thetaXW (degrees): angle of the track projected onto the (wire, drift) plane.
// wire_pitch is the wire spacing (0.3 cm for all three ICARUS planes).
inline double ThetaXW(double dirx, double pitch, double wire_pitch = 0.3)
{
    return std::atan(dirx * pitch / wire_pitch) * 180. / TMath::Pi();
}

// Returns the ratio_integral scale factor for a given hit.
// x        : hit x position in detector coordinates [cm]
// dirx     : track direction x component (from trk.dir.x)
// pitch    : effective hit pitch [cm] (from calo.points[i].pitch)
// tpc      : TPC index (0-3)
// plane    : wire plane (0-2; use 2 for collection)
// run_period: 2 or 4
inline double GetWireModScaleFactor(double x, double dirx, double pitch,
                                     int tpc, int plane, int run_period)
{
    return wiremod_table.scale_factor(ThetaXW(dirx, pitch), x, tpc, plane, run_period);
}

// --------------------------------------------------------------------------

inline Var MakeMuonChi2MuCorrected(double f)
{
    return Var([f](const caf::SRSliceProxy* islc) -> double {
        int ipfp_mu = kIcarus202401MuonIdx(islc);
        if (ipfp_mu < 0) return -1.;
        const auto& trk = islc->reco.pfp[ipfp_mu].trk;
        double dirx = trk.dir.x;
        auto chi2 = chi2_corrected_calculator.calculate_chi2_corrected(trk.calo[2], dirx, f);
        return chi2.chi2_muon;
    });
}

inline Var MakeLeadingProtonChi2ProtonCorrected(double f)
{
    return Var([f](const caf::SRSliceProxy* islc) -> double {
        const auto ps = kIcarus202401RecoProtonP(islc);
        if (!ps.size()) return -1.;
        double p = *std::max_element(ps.begin(), ps.end());
        for (const auto& pfp : islc->reco.pfp) {
            if (pfp.trk.rangeP.p_proton != p) continue;
            double dirx = pfp.trk.dir.x;
            return chi2_corrected_calculator
                .calculate_chi2_corrected(pfp.trk.calo[2], dirx, f)
                .chi2_proton;
        }
        return -1.;
    });
}

} // namespace dedxsyst
} // namespace ana
