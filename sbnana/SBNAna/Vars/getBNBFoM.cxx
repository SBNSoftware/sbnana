//! ////////////////////////////////////////////////////////////////////////////
//! @file: getBNBFoM.cxx
//! @author: Jacob Smith (smithja)
//! @email: jacob.a.smith@stonybrook.edu
//!
//! Special thanks to Zarko Pavlovic (zarko) and Joseph Zennamo (jaz8600) for
//! helping me understand the MicroBooNE code which this script is based off of.
//!
//! @note This code is almost entirely based off MicroBooNE's getFOM2 code in 
//! the ubraw GitHub repo. How we calculate our Figure of Merit (FoM) is almost
//! identical to what's done in getFOM2.cxx. See below for further details.
//!
//! @details Calculates the Booster Neutrino Beam (BNB) Figure of Merit (FoM),
//! which is a score on the interval [0, 1] that acts as a measure of the 
//! geometric overlap of the BNB with the beamline's nuclear target. There are 
//! separate cuts placed on the protons per pulse (PPP) and horn current in the 
//! BNBVars<.h, .cxx> and BNBQualityCuts<.h, .cxx> files. The PPP––via 
//! kSpillTOR875 (and secondarily kSpillTOR876)-––is the only quantity from the
//! suite referenced in the above files that makes its way into FoM calculations.
//! ////////////////////////////////////////////////////////////////////////////
#include "getBNBFoM.h"

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <math.h>
#include <stdexcept>

#include "TROOT.h"
#include "TH1D.h"
#include "TF1.h"


/** @struct OffsetData Holds the timestamped calibration positions of Beam
 * Position Monitors.
 */
struct OffsetData {
    std::vector<int> Times; //! non-leap seconds since UNIX epoch
    std::vector<Double_t> Vals; //! transverse distance from BNB beamline [mm]
};

std::unordered_map<std::string, OffsetData> offsets = {
    {"HP875", {
        { 1420092000, 1574190445, 1576014670, 1588603117, 1605634657, 
            1606510357, 1608593857, 1668802510, 1711030691 },
        { -3.40, -1.20, -1.20, -3.33, -1.20, 0.75, -2.40, -5.53, -3.42 }
    }},
    {"VP875", {
        { 1420092000, 1574190445, 1576014670, 1588603117, 1605634657, 
            1606510357, 1608070057, 1608593857, 1612921957 },
        { 1.48, 1.30, 1.30, 1.43, 1.30, -0.34, -3.60, 2.40, -0.70 }
    }},
    {"HPTG1", {
        { 1420092000, 1576014670, 1606510357, 1677887110, 1711030691 },
        { 0.46, 0.0, -2.15, -1.16, -0.03 }
    }},
    {"VPTG1", {
        { 1420092000, 1576014670, 1588603117, 1605706657, 1606510357,
            1608070057, 1608593857 },
        { 0.39, -0.10, 0.0, -0.10, -4.53, -7.50, -4.70 }
    }},
    {"HPTG2", {
        { 1420092000, 1574190444, 1576014670, 1588603117, 1605706657, 
            1606510357, 1608593857, 1668802510 },
        { 1.07, -0.40, -0.40, 0.0, -0.40, -2.31, 0.0, -0.84 }
    }},
    {"VPTG2", {
        { 1420092000, 1572899153, 1574190445, 1576014670, 1588603117,
            1605706657, 1606510357, 1612798657, 1612921957 },
        { 1.84, -3.15, 0.10, 0.10, 0.0, 0.10, 0.82, 0.0, 0.80 }
    }},
    {"VP873", {
        { 1615250280, 1734021000, 1734034620 },
        { 3.4, 12, 2.8 }
    }}
};

//! Position of beam-monitoring device along the direction of the 
//! beamline, i.e. z-position [mm].
Double_t VP873ZPos         = 1911.53656;
Double_t mw875ZPos         = 2018.20648;
Double_t HP875ZPos         = 2021.16104;
Double_t VP875ZPos         = 2023.193205;
Double_t mw876ZPos         = 2029.66766;
Double_t HPTG1ZPos         = 2048.33267;
Double_t HPTG2ZPos         = 2052.40662;
Double_t targetCenterZPos  = 2068.70895;

//! Numbers for beam width measurements from M875BB and M876BB.
//! These are used in 0th, 1st, and 2nd order expansions when
//! determining beam width quantities.
//! Study done by Zarko Pavlovic (zarko) to obtain these numbers.
//! How he or others arived at these numbers and what they represent is unknown.
//! But given the low 2nd-order expansion coefficents, it's unlikely this study
//! will ever need to be repeated.
std::vector<Double_t> p875x = { 0.431857, 0.158077, 0.00303551};
std::vector<Double_t> p875y = { 0.279128, 0.337048, 0};
std::vector<Double_t> p876x = { 0.166172, 0.30999,  -0.00630299};
std::vector<Double_t> p876y = { 0.13425,  0.580862, 0};

/** @fn expandWidth()
 * @brief Utilizing BNB studies (see notes around p876x, etc.), calculate
 * the BNB width at the nuclear target up to an N-th order expansion
 * @return BNB width at the nuclear target
 */
static std::optional<Double_t> expandWidth(const Double_t width, const std::vector<Double_t>& coeff, const size_t N) {
    if (N > coeff.size()) {
        std::cerr << "[WARN] Requested expansion order " << N 
                  << " exceeds number of coefficients (" << coeff.size() 
                  << ")." << std::endl;
        return std::nullopt;
    }
    Double_t result = coeff[0];
    for (size_t i = 1; i < N; ++i) {
        result += coeff[i] * std::pow(width, i);
    }
    return result;
}

/** @fn isReasonable()
 * @brief tests if a BNB width is within the given bounds
 * @return boolean indicating if given width is within given bounds
 */
static bool isReasonable( const Double_t width, const Double_t lower, const Double_t upper) {
    if ( width >= lower && width <= upper) return true;
    return false;
}

/** @fn getValidBPMCalibVal()
 * @brief Return the latest calibration value before or at a given timestamp for
 * a specific device.
 * 
 * @param timestamp UNIX-formatted timestamp [s].
 * @param device The name of the beam monitoring device.
 * 
 * @return Optional value: for given timestamp, most recent valid calibration 
 * value [mm] (or std::nullopt).
 */
std::optional<Double_t> getValidBPMCalibVal(const int& timestamp,
                                        const std::string& device) {
    const auto& times = offsets.at(device).Times;
    const auto& vals  = offsets.at(device).Vals;

    //! Find the first index where times[i] > timestamp
    //! The last time <= timestamp is then at index (i - 1), if i != 0
    auto it = std::upper_bound(times.begin(), times.end(), timestamp);

    if (it == times.begin()) {
        return std::nullopt; //! No value at or before timestamp
    } 

    int idx = std::distance(times.begin(), it) - 1;
    return vals[idx]; //! Latest valid value
}

/** @fn getValidBPMCalibTime()
 * @brief Return the latest timestamp of a calibration value before or at a 
 * given timestamp for a specific device.
 * 
 * @param timestamp UNIX-formatted timestamp [s].
 * @param device The name of the beam monitoring device.
 * 
 * @return Optional value: for given timestamp, most recent timestamp of a 
 * valid calibration value [mm] (or std::nullopt).
 */
std::optional<Double_t> getValidBPMCalibTime(const int& timestamp,
                                        const std::string& device) {
    const auto& times = offsets.at(device).Times;

    //! Find the first index where times[i] > timestamp
    //! The last time <= timestamp is then at index (i - 1), if i != 0
    auto it = std::upper_bound(times.begin(), times.end(), timestamp);

    if (it == times.begin()) {
        return std::nullopt; //! No valid value at or before timestamp
    } 

    int idx = std::distance(times.begin(), it) - 1;
    return times[idx]; //! Latest timestamp for a valid value
}

/** @fn swimBNB()
 * @brief Propagate beam centroid and sigma matrix through transfer matrices.
 *
 * @param centroid1 Input centroid (6-vector)
 * @param sigma1 Input sigma matrix (6x6)
 * @param xferc Transfer matrix for centroid
 * @param xfers Transfer matrix for sigma
 * @param cx Output centroid x [mm]
 * @param cy Output centroid y [mm]
 * @param sx Output beam sigma-x [mm]
 * @param sy Output beam sigma-y [mm]
 * @param rho Output beam correlation coefficient
 * 
 * @post Beam centriod and sigma parameters as well as beam correlation 
 * coefficent are modifed in the body of calcFoM() or @calcFoM2()
 */
void swimBNB(const Double_t centroid1[6], 
                const Double_t sigma1[6][6],
                const Double_t xferc[6][6], 
                const Double_t xfers[6][6],
                Double_t& cx, 
                Double_t& cy, 
                Double_t& sx, 
                Double_t& sy, 
                Double_t& rho) {
    Double_t centroid2[6] = {0};
    for (unsigned int i = 0; i < 6; ++i)
        for (unsigned int j = 0; j < 6; ++j)
            centroid2[i] += xferc[i][j] * centroid1[j];

    cx = centroid2[0];
    cy = centroid2[2];

    Double_t sigma2[6][6] = {{0}};
    for (unsigned int i = 0; i < 6; ++i)
        for (unsigned int j = 0; j < 6; ++j)
            for (unsigned int k = 0; k < 6; ++k)
                for (unsigned int m = 0; m < 6; ++m)
                    sigma2[i][m] += xfers[i][j] * sigma1[j][k] * xfers[m][k];

    sx = std::sqrt(sigma2[0][0]) * 1000.0;
    sy = std::sqrt(sigma2[2][2]) * 1000.0;
    rho = sigma2[0][2] / std::sqrt(sigma2[0][0] * sigma2[2][2]);
}

/** @fn funcIntBivar()
 * @brief Compute overlap integral of a 2D bivariate Gaussian with a cylinder.
 * 
 * @param cx Beam centroid x-position [mm]
 * @param cy Beam centroid y-position [mm]
 * @param sx Beam sigma-x [mm]
 * @param sy Beam sigma-y [mm]
 * @param rho Beam correlation coefficient
 * 
 * @return log10(1 - overlap fraction), or -10000 if integral is unphysical
 */
Double_t funcIntBivar(const Double_t cx, 
                        const Double_t cy,
                        const Double_t sx, 
                        const Double_t sy,
                        const Double_t rho) {
    const Double_t dbin = 0.1, r = 4.75, rr = r * r;
    const Double_t rho2 = rho * rho;
    const int imax = static_cast<int>(round((2.0 * r) / dbin));
    const int jmax = static_cast<int>(round((2.0 * r) / dbin));
    Double_t sum = 0.0;

    for (unsigned int i = 0; i <= imax; ++i) {
        Double_t x = -r + i * dbin;
        for (unsigned int j = 0; j <= jmax; ++j) {
            Double_t y = -r + j * dbin;
            if (x * x + y * y < rr) {
                Double_t tx = (x + cx) / sx;
                Double_t ty = (y + cy) / sy;
                Double_t z = tx * tx - 2.0 * rho * tx * ty + ty * ty;
                sum += std::exp(-z / (2.0 * (1.0 - rho2)));
            }
        }
    }

    sum *= dbin * dbin / (2.0 * 3.14159 * sx * sy * std::sqrt(1.0 - rho2));
    return (sum >= 1.0) ? -10000. : std::log10(1.0 - sum);
}

/** @fn calcFoM()
 * @brief Calculate the Figure of Merit for BNB alignment WITHOUT use of 
 * multi-wire data.
 * 
 * @details This function uses the beam's centroid and sigma matrix to propagate
 * ("swim") the beam envelope to different locations (upstream, center, 
 * downstream) of the target, computes overlap integrals with a cylindrical 
 * target, and forms a weighted sum.
 * 
 * @note calcFoM() is the same as calcFoM2() but with scalex = scaley = 1 
 * 
 * @param horPos BNB's horizontal position at nuclear target [mm]
 * @param horAng BNB's horizontal angle at nuclear target [rad]
 * @param verPos BNB's vertical position at nuclear target [mm]
 * @param verAng BNB's vertical angle at nuclear target [rad]
 * @param PPP Protons per pulse (used for emittance and momentum spread)
 * 
 * @return Double representing log10(1 - overlap-fraction), or -10000 if invalid.
 */
Double_t calcFoM(const Double_t horPos, 
                    const Double_t horAng, 
                    const Double_t verPos, 
                    const Double_t verAng,
                    const Double_t PPP) {
    //! Beam Twiss parameters and dispersions from MiniBooNE AnalysisFramework
    //! Code from DQ_BeamLine_twiss_init.F
    const Double_t bx = 4.68, ax = 0.0389, gx = (1 + ax * ax) / bx;
    const Double_t nx = 0.0958, npx = -0.0286;
    const Double_t by = 59.12, ay = 2.4159, gy = (1 + ay * ay) / by;
    const Double_t ny = 0.4577, npy = -0.0271;

    //! LIKELY emittance and momentum spread as functions of protons-per-pulse
    //! Code from DQ_BeamLine_make_tgt_fom2.F
    const Double_t ex = 0.1775E-06 + 0.1827E-07 * PPP;
    const Double_t ey = 0.1382E-06 + 0.2608E-08 * PPP;
    const Double_t dp = 0.4485E-03 + 0.6100E-04 * PPP;

    Double_t tex = ex, tey = ey, tdp = dp;

    //! Code from DQ_BeamLine_beam_init.F
    Double_t sigma1[6][6] = {{0}};
    Double_t centroid1[6] = {horPos, horAng, verPos, verAng, 0.0, 0.0};

    sigma1[5][5] = tdp * tdp;
    sigma1[0][0] = tex * bx + nx * nx * tdp * tdp;
    sigma1[0][1] = -tex * ax + nx * npx * tdp * tdp;
    sigma1[1][1] = tex * gx + npx * npx * tdp * tdp;
    sigma1[0][5] = nx * tdp * tdp;
    sigma1[1][5] = npx * tdp * tdp;
    sigma1[1][0] = sigma1[0][1];
    sigma1[5][0] = sigma1[0][5];
    sigma1[5][1] = sigma1[1][5];

    sigma1[2][2] = tey * by + ny * ny * tdp * tdp;
    sigma1[2][3] = -tey * ay + ny * npy * tdp * tdp;
    sigma1[3][3] = tey * gy + npy * npy * tdp * tdp;
    sigma1[2][5] = ny * tdp * tdp;
    sigma1[3][5] = npy * tdp * tdp;
    sigma1[3][2] = sigma1[2][3];
    sigma1[5][2] = sigma1[2][5];
    sigma1[5][3] = sigma1[3][5];

    const Double_t begtocnt[6][6] = {
        {0.65954, 0.43311, 0.00321, 0.10786, 0.00000, 1.97230},
        {0.13047, 1.60192, 0.00034, 0.00512, 0.00000, 1.96723},
        {-0.00287, -0.03677, -0.35277, -4.68056, 0.00000, 0.68525},
        {-0.00089, -0.00430, -0.17722, -5.18616, 0.00000, 0.32300},
        {-0.00104, 0.00232, -0.00001, -0.00224, 1.00000, -0.00450},
        {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 1.00000}
    };

    Double_t cnttoups[6][6] = {{0}};
    Double_t cnttodns[6][6] = {{0}};
    Double_t identity[6][6] = {{0}};
    Double_t begtoups[6][6] = {{0}};
    Double_t begtodns[6][6] = {{0}};

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j) {
            identity[i][j] = (i == j);
            cnttoups[i][j] = (i == j);
            cnttodns[i][j] = (i == j);
        }
    }
    cnttoups[0][1] = cnttoups[2][3] = -0.35710;
    cnttodns[0][1] = cnttodns[2][3] = +0.35710;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j) {
            for (unsigned int k = 0; k < 6; ++k) {
                begtoups[i][k] += cnttoups[i][j] * begtocnt[j][k];
                begtodns[i][k] += cnttodns[i][j] * begtocnt[j][k];
            }
        }
    }

    Double_t cx, cy, sx, sy, rho;

    //! swim to upstream of target
    swimBNB(centroid1, sigma1, cnttoups, begtoups, cx, cy, sx, sy, rho);
    Double_t fom_a = funcIntBivar(cx, cy, sx, sy, rho);

    //! swim to center of target
    swimBNB(centroid1, sigma1, identity, begtocnt, cx, cy, sx, sy, rho);
    Double_t fom_b = funcIntBivar(cx, cy, sx, sy, rho);

    //! swim to downstream of target
    swimBNB(centroid1, sigma1, cnttodns, begtodns, cx, cy, sx, sy, rho);
    Double_t fom_c = funcIntBivar(cx, cy, sx, sy, rho);

    if (fom_a <= -10000. 
        || fom_b <= -10000. 
        || fom_c <= -10000.) {
            return -10000.;
    }

    return fom_a * 0.6347 + fom_b * 0.2812 + fom_c * 0.0841;
}

/** @fn calcFoM2()
 * @brief Calculate the Figure of Merit for BNB alignment WITH the use of 
 * multi-wire data.
 * 
 * @details This function uses the beam's centroid and sigma matrix to propagate
 * ("swim") the beam envelope to different locations (upstream, center, 
 * downstream) of the target, computes overlap integrals with a cylindrical 
 * target, and forms a weighted sum.
 * 
 * @param horPos BNB's horizontal position at nuclear target [mm]
 * @param horAng BNB's horizontal angle at nuclear target [rad]
 * @param verPos BNB's vertical position at nuclear target [mm]
 * @param verAng BNB's vertical angle at nuclear target [rad]
 * @param PPP Protons per pulse (used for emittance and momentum spread)
 * @param tgtHorWidth Horizontal width of the BNB at the nuclear target [mm]
 * @param tgtVerWidth Vertical width of the BNB at the nuclear target [mm]
 * 
 * @return Double representing log10(1 - overlap-fraction), or -10000 if invalid.
 */
Double_t calcFoM2(const Double_t horPos, 
                    const Double_t horAng, 
                    const Double_t verPos, 
                    const Double_t verAng,
                    const Double_t PPP, 
                    const Double_t tgtHorWidth, 
                    const Double_t tgtVerWidth) {
    //! Beam Twiss parameters and dispersions from MiniBooNE AnalysisFramework
    //! Code from DQ_BeamLine_twiss_init.F
    const Double_t bx = 4.68, ax = 0.0389, gx = (1 + ax * ax) / bx;
    const Double_t nx = 0.0958, npx = -0.0286;
    const Double_t by = 59.12, ay = 2.4159, gy = (1 + ay * ay) / by;
    const Double_t ny = 0.4577, npy = -0.0271;

    //! LIKELY emittance and momentum spread as functions of protons-per-pulse
    //! Code from DQ_BeamLine_make_tgt_fom2.F
    const Double_t ex = 0.1775E-06 + 0.1827E-07 * PPP;
    const Double_t ey = 0.1382E-06 + 0.2608E-08 * PPP;
    const Double_t dp = 0.4485E-03 + 0.6100E-04 * PPP;

    Double_t tex = ex, tey = ey, tdp = dp;

    //! Code from DQ_BeamLine_beam_init.F
    Double_t sigma1[6][6] = {{0}};
    Double_t centroid1[6] = {horPos, horAng, verPos, verAng, 0.0, 0.0};

    sigma1[5][5] = tdp * tdp;
    sigma1[0][0] = tex * bx + nx * nx * tdp * tdp;
    sigma1[0][1] = -tex * ax + nx * npx * tdp * tdp;
    sigma1[1][1] = tex * gx + npx * npx * tdp * tdp;
    sigma1[0][5] = nx * tdp * tdp;
    sigma1[1][5] = npx * tdp * tdp;
    sigma1[1][0] = sigma1[0][1];
    sigma1[5][0] = sigma1[0][5];
    sigma1[5][1] = sigma1[1][5];

    sigma1[2][2] = tey * by + ny * ny * tdp * tdp;
    sigma1[2][3] = -tey * ay + ny * npy * tdp * tdp;
    sigma1[3][3] = tey * gy + npy * npy * tdp * tdp;
    sigma1[2][5] = ny * tdp * tdp;
    sigma1[3][5] = npy * tdp * tdp;
    sigma1[3][2] = sigma1[2][3];
    sigma1[5][2] = sigma1[2][5];
    sigma1[5][3] = sigma1[3][5];

    const Double_t begtocnt[6][6] = {
        {0.65954, 0.43311, 0.00321, 0.10786, 0.00000, 1.97230},
        {0.13047, 1.60192, 0.00034, 0.00512, 0.00000, 1.96723},
        {-0.00287, -0.03677, -0.35277, -4.68056, 0.00000, 0.68525},
        {-0.00089, -0.00430, -0.17722, -5.18616, 0.00000, 0.32300},
        {-0.00104, 0.00232, -0.00001, -0.00224, 1.00000, -0.00450},
        {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 1.00000}
    };

    Double_t cnttoups[6][6] = {{0}};
    Double_t cnttodns[6][6] = {{0}};
    Double_t identity[6][6] = {{0}};
    Double_t begtoups[6][6] = {{0}};
    Double_t begtodns[6][6] = {{0}};

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j) {
            identity[i][j] = (i == j);
            cnttoups[i][j] = (i == j);
            cnttodns[i][j] = (i == j);
        }
    }
    cnttoups[0][1] = cnttoups[2][3] = -0.35710;
    cnttodns[0][1] = cnttodns[2][3] = +0.35710;

    for (unsigned int i = 0; i < 6; ++i) {
        for (unsigned int j = 0; j < 6; ++j) {
            for (unsigned int k = 0; k < 6; ++k) {
                begtoups[i][k] += cnttoups[i][j] * begtocnt[j][k];
                begtodns[i][k] += cnttodns[i][j] * begtocnt[j][k];
            }
        }
    }

    Double_t cx, cy, sx, scalex, sy, scaley, rho;

    //! swim to upstream of target
    swimBNB(centroid1, sigma1, cnttoups, begtoups, cx, cy, sx, sy, rho);
    scalex = tgtHorWidth / sx;
    scaley = tgtVerWidth / sy;
    Double_t fom_a = funcIntBivar(cx, cy, sx*scalex, sy*scaley, rho);

    //! swim to center of target
    swimBNB(centroid1, sigma1, identity, begtocnt, cx, cy, sx, sy, rho);
    scalex = tgtHorWidth / sx;
    scaley = tgtVerWidth / sy;
    Double_t fom_b = funcIntBivar(cx, cy, sx*scalex, sy*scaley, rho);

    //! swim to downstream of target
    swimBNB(centroid1, sigma1, cnttodns, begtodns, cx, cy, sx, sy, rho);
    scalex = tgtHorWidth / sx;
    scaley = tgtVerWidth / sy;
    Double_t fom_c = funcIntBivar(cx, cy, sx*scalex, sy*scaley, rho);

    if (fom_a <= -10000. 
        || fom_b <= -10000. 
        || fom_c <= -10000.) {
            return -10000.;
    }

    return fom_a * 0.6347 + fom_b * 0.2812 + fom_c * 0.0841;
}

/** @fn getBNBFoM()
 * @brief Identifies and calculates necessary quantities to calculate a Figure
 * of Merit then returns the Figure of Merit.
 * 
 * @param spillTimeSec readout timestamp for given spill [non-leap seconds
 * since UNIX epoch (Jan 1, 1970 at midnight UTC/GMT)]
 * @param TOR Protons per pulse (used for emittance and momentum spread)
 * @param HP875 Horizontal position of BNB upstream of nuclear target [mm]
 * @param HPTG1 Horizontal position of BNB upstream of nuclear target [mm]
 * @param HPTG2 Horizontal position of BNB upstream of nuclear target [mm]
 * @param VP873 Vertical position of BNB upstream of nuclear target [mm]
 * @param VP875 Vertical position of BNB upstream of nuclear target [mm]
 * 
 * Order of horizontal BPMs: HP875, HPTG1, HPTG2
 * Order of vertical BPMs: VP873, VP875
 * 
 * @return BNB Figure of Merit not utilizing beam width data from multi-wire
 * devices
 */
Double_t getBNBFoM( const Double_t spillTimeSec,
                      const Double_t TOR,
                      const Double_t HP875, 
                      const Double_t HPTG1,
                      const Double_t HPTG2,
                      const Double_t VP873,
                      const Double_t VP875) {
    /** 
     * HPTG1 and HPTG2 calibration positions are expected to be complementary, 
     * i.e. when one has a valid calibration position the other does not. To
     * determine which of these BPMs to use, we pick the one that has the
     * most recently updated (and valid) calibration position. Furthermore, we
     * prioritize the use of HPTG2 if there is ever a spillTimeSec where
     * HPTG1TimeDiff == HPTG2TimeDiff since HPTG2 is marginally closer to the
     * nuclear target. This helps to give the most accurate horAng and 
     * tgtHorPos when projecting from HP875.
     */
    Double_t HP875TimeDiff = -1.;
    Double_t HPTG1TimeDiff = -1.;
    Double_t HPTG2TimeDiff = -1.;
    Double_t VP873TimeDiff = -1.;
    Double_t VP875TimeDiff = -1.;

    std::optional<Double_t> HP875CalibTime = getValidBPMCalibTime( spillTimeSec, "HP875");
    std::optional<Double_t> HPTG1CalibTime = getValidBPMCalibTime( spillTimeSec, "HPTG1");
    std::optional<Double_t> HPTG2CalibTime = getValidBPMCalibTime( spillTimeSec, "HPTG2");

    std::optional<Double_t> VP873CalibTime = getValidBPMCalibTime( spillTimeSec, "VP873");
    std::optional<Double_t> VP875CalibTime = getValidBPMCalibTime( spillTimeSec, "VP875");

    if (HP875CalibTime.has_value()) {
        HP875TimeDiff = spillTimeSec - HP875CalibTime.value();
    }

    if (HPTG1CalibTime.has_value()) {
        HPTG1TimeDiff = spillTimeSec - HPTG1CalibTime.value();
    }

    if (HPTG2CalibTime.has_value()) {
        HPTG2TimeDiff = spillTimeSec - HPTG2CalibTime.value();
    }

    if (VP873CalibTime.has_value()) {
        VP873TimeDiff = spillTimeSec - VP873CalibTime.value();
    }

    if (VP875CalibTime.has_value()) {
        VP875TimeDiff = spillTimeSec - VP875CalibTime.value();
    }


    /** @warning
     * If both HPTG1TimeDiff and HPTG2TimeDiff are negative (very unlikely),
     * the earliest known HPTG1 and HPTG2 calibration positions were logged
     * after the current spillTimeSec. This prevents us from calculating a
     * Figure of Merit since we only have one horizontal BPM and cannot
     * calculate horAng and tgtHorPos. Similarly, if we have at least one of
     * VP873TimeDiff and VP875TimeDiff less than zero, we'll only have one
     * vertical BPM and will not be able to calculate a Figure of Merit. Since
     * HP875 is theoretically valid for all ICARUS timestamps, we also assert
     * we could not calculate a (meaningful) Figure of Merit if HP875 *ever*
     * has a negative time difference for a given spillTimeSec timestamp.
     */
    if ( (HPTG1TimeDiff < 0 && HPTG2TimeDiff < 0)
            || (HP875TimeDiff < 0 || VP873TimeDiff < 0 || VP875TimeDiff < 0)) { 
        return -2.;
    }
    
    //! Get most recent calibration position of Beam Position Monitor (BPM) 
    //! devices and calculate offset from calibration position.
    std::optional<Double_t> HP875CalibVal = getValidBPMCalibVal( spillTimeSec, "HP875");
    Double_t HP875Offset = HP875 - HP875CalibVal.value();

    Double_t HPTGXOffset = -999; //! 'X' is stand-in 
    Double_t HPTGXZPos = -999;   //! char for '1' or '2'
    if (HPTG2TimeDiff <= HPTG1TimeDiff) { //! prefer HPTG2 if time diffs match
        std::optional<Double_t> HPTG2CalibVal = getValidBPMCalibVal( spillTimeSec, "HPTG2");
        HPTGXOffset = HPTG2 - HPTG2CalibVal.value();        
        HPTGXZPos = HPTG2ZPos;
    }
    else {
        std::optional<Double_t> HPTG1CalibVal = getValidBPMCalibVal( spillTimeSec, "HPTG1");
        HPTGXOffset = HPTG1 - HPTG1CalibVal.value();        
        HPTGXZPos = HPTG1ZPos;
    }

    std::optional<Double_t> VP873CalibVal = getValidBPMCalibVal( spillTimeSec, "VP873");
    Double_t VP873Offset = VP873 - VP873CalibVal.value();    
    
    std::optional<Double_t> VP875CalibVal = getValidBPMCalibVal( spillTimeSec, "VP875");
    Double_t VP875Offset = VP875 - VP875CalibVal.value();
    /**
     * Calculate angles between horizontal/vertical BPMs and project the BNB 
     * to the transverse plane intersecting the center of the nuclear target 
     * (along the beamline direction, i.e. z). We project the horizontal
     * position from HP875 to the nuclear target since HPTG1 and HPTG2 are close 
     * to the front of the target. We project the vertical position from VP875
     * since both VP873 and VP875 are sufficiently far from the nuclear target
     * and VP875 is closer to the nuclear target.
     * 
     * @note angle calculations use the small angle approximation
     * 
     * Order of horizontal BPMs: HP875, HPTG1, HPTG2
     * Order of vertical BPMs: VP873, VP875
     */
    Double_t horAng = (HPTGXOffset - HP875Offset) / (HPTGXZPos - HP875ZPos);
    horAng = atan(horAng);
    Double_t tgtHorPos = HP875Offset + horAng * (targetCenterZPos - HP875ZPos);

    Double_t verAng = (VP875Offset - VP873Offset) / (VP875ZPos - VP873ZPos);
    verAng = atan(verAng);
    Double_t tgtVerPos = VP873Offset + horAng * (targetCenterZPos - HP875ZPos);

    /**
     * Calculating a Figure of Merit without the use of multi-wire data amounts
     * to setting scalex and scaley both equal to 1 (see documentation on
     * getBNBFoM() and getBNBFoM2() for more details).
     * scalex = tgtHorWidth / sx, i.e. the x-scale is the width of the BNB at 
     * the center of the nuclear target divided by the width of the BNB at at 
     * one of the positions you can "swim" to (see swimBNB for details on this).
     * In general, we have that tgtHorWidth = sx for a given timestamp. At this
     * timestamp, we can "swim" to different locations to get different sx, but
     * tgtHorWidth is constant at a given timestamp. We set tgtHorWidth = sx = 1 
     * for this Figure of Merit calculation that does not make use of the 
     * multi-wire data. tgtHorWidth = sx = 1 is done inside calcFoM(), and thus 
     * we do not pass BNB width data to the calcFoM() below. For a Figure of 
     * Merit that makes use of BNB width data, please see getBNBFoM2() and
     * calcFoM2().
     *
     * @note One could do studies to come up with a better value for 
     * tgtHorWidth = sx by, e.g., visually inspecting the (hopefully) gaussian 
     * profiles of the multi-wire devices. However, since the goal in SBN is to 
     * use a Figure of Merit that *does* use multi-wire data, we have not 
     * concerned ourselves too much with what constant value we set 
     * tgtHorWidth = sx.
     */
    return calcFoM( tgtHorPos, horAng, tgtVerPos, verAng, TOR);
}

/** @fn getBNBFoM2()
 * @brief Identifies and calculates necessary quantities to calculate a Figure
 * of Merit then returns the Figure of Merit.
 * 
 * @note If there is not good multi-wire fit parameters, this function defaults
 * to returning a Figure of Merit that does not use multi-wire data, i.e.
 * calcFoM().
 * 
 * @param spillTimeSec readout timestamp for given spill [non-leap seconds
 * since UNIX epoch (Jan 1, 1970 at midnight UTC/GMT)]
 * @param TOR Protons per pulse (used for emittance and momentum spread)
 * @param HP875 Horizontal position of BNB upstream of nuclear target [mm]
 * @param HPTG1 Horizontal position of BNB upstream of nuclear target [mm]
 * @param HPTG2 Horizontal position of BNB upstream of nuclear target [mm]
 * @param VP873 Vertical position of BNB upstream of nuclear target [mm]
 * @param VP875 Vertical position of BNB upstream of nuclear target [mm]
 * @param MW875HorWidth Horizontal width of the BNB upstream of target [mm]
 * @param MW875VerWidth Vertical width of the BNB upstream of target [mm]
 * @param MW876HorWidth Horizontal width of the BNB upstream of target [mm]
 * @param MW876VerWidth Vertical width of the BNB upstream of target [mm]
 * 
 * Order of horizontal BPMs: HP875, HPTG1, HPTG2
 * Order of vertical BPMs: VP873, VP875
 * Order of multi-wire readout (MWR) devices: MW875, MW876
 * 
 * @return BNB Figure of Merit utilizing beam width data from multi-wire
 * devices.
 */
Double_t getBNBFoM2( const Double_t spillTimeSec,
                      const Double_t TOR,
                      const Double_t HP875, 
                      const Double_t HPTG1,
                      const Double_t HPTG2,
                      const Double_t VP873,
                      const Double_t VP875,
                      const Double_t MW875HorWidth,
                      const Double_t MW875VerWidth,
                      const Double_t MW876HorWidth,
                      const Double_t MW876VerWidth) {
    /** 
     * HPTG1 and HPTG2 calibration positions are expected to be complementary, 
     * i.e. when one has a valid calibration position the other does not. To
     * determine which of these BPMs to use, we pick the one that has the
     * most recently updated (and valid) calibration position. Furthermore, we
     * prioritize the use of HPTG2 if there is ever a spillTimeSec where
     * HPTG1TimeDiff == HPTG2TimeDiff since HPTG2 is marginally closer to the
     * nuclear target. This helps to give the most accurate horAng and 
     * tgtHorPos when projecting from HP875.
     */
    Double_t HP875TimeDiff = -1.;
    Double_t HPTG1TimeDiff = -1.;
    Double_t HPTG2TimeDiff = -1.;
    Double_t VP873TimeDiff = -1.;
    Double_t VP875TimeDiff = -1.;

    std::optional<Double_t> HP875CalibTime = getValidBPMCalibTime( spillTimeSec, "HP875");
    std::optional<Double_t> HPTG1CalibTime = getValidBPMCalibTime( spillTimeSec, "HPTG1");
    std::optional<Double_t> HPTG2CalibTime = getValidBPMCalibTime( spillTimeSec, "HPTG2");

    std::optional<Double_t> VP873CalibTime = getValidBPMCalibTime( spillTimeSec, "VP873");
    std::optional<Double_t> VP875CalibTime = getValidBPMCalibTime( spillTimeSec, "VP875");

    if (HP875CalibTime.has_value()) {
        HP875TimeDiff = spillTimeSec - HP875CalibTime.value();
    }

    if (HPTG1CalibTime.has_value()) {
        HPTG1TimeDiff = spillTimeSec - HPTG1CalibTime.value();
    }

    if (HPTG2CalibTime.has_value()) {
        HPTG2TimeDiff = spillTimeSec - HPTG2CalibTime.value();
    }

    if (VP873CalibTime.has_value()) {
        VP873TimeDiff = spillTimeSec - VP873CalibTime.value();
    }

    if (VP875CalibTime.has_value()) {
        VP875TimeDiff = spillTimeSec - VP875CalibTime.value();
    }

    /** @warning
     * If both HPTG1TimeDiff and HPTG2TimeDiff are negative (very unlikely),
     * the earliest known HPTG1 and HPTG2 calibration positions were logged
     * after the current spillTimeSec. This prevents us from calculating a
     * Figure of Merit since we only have one horizontal BPM and cannot
     * calculate horAng and tgtHorPos. Similarly, if we have at least one of
     * VP873TimeDiff and VP875TimeDiff less than zero, we'll only have one
     * vertical BPM and will not be able to calculate a Figure of Merit. Since
     * HP875 is theoretically valid for all ICARUS timestamps, we also assert
     * we could not calculate a (meaningful) Figure of Merit if HP875 *ever*
     * has a negative time difference for a given spillTimeSec timestamp.
     */
    if ( (HPTG1TimeDiff < 0 && HPTG2TimeDiff < 0)
            || (HP875TimeDiff < 0 || VP873TimeDiff < 0 || VP875TimeDiff < 0)) { 
        return -2.;
    }
    
    //! Get most recent calibration position of Beam Position Monitor (BPM) 
    //! devices and calculate offset from calibration position.
    std::optional<Double_t> HP875CalibVal = getValidBPMCalibVal( spillTimeSec, "HP875");
    Double_t HP875Offset = HP875 - HP875CalibVal.value();

    Double_t HPTGXOffset = -999; //! 'X' is stand-in 
    Double_t HPTGXZPos = -999;   //! char for '1' or '2'
    if (HPTG2TimeDiff <= HPTG1TimeDiff) { //! prefer HPTG2 if time diffs match
        std::optional<Double_t> HPTG2CalibVal = getValidBPMCalibVal( spillTimeSec, "HPTG2");
        HPTGXOffset = HPTG2 - HPTG2CalibVal.value();  
        HPTGXZPos = HPTG2ZPos;
    }
    else {
        std::optional<Double_t> HPTG1CalibVal = getValidBPMCalibVal( spillTimeSec, "HPTG1");
        HPTGXOffset = HPTG1 - HPTG1CalibVal.value();
        HPTGXZPos = HPTG1ZPos;
    }

    std::optional<Double_t> VP873CalibVal = getValidBPMCalibVal( spillTimeSec, "VP873");
    Double_t VP873Offset = VP873 - VP873CalibVal.value();
    
    std::optional<Double_t> VP875CalibVal = getValidBPMCalibVal( spillTimeSec, "VP875");
    Double_t VP875Offset = VP875 - VP875CalibVal.value();
    /**
     * Calculate angles between horizontal/vertical BPMs and project the BNB 
     * to the transverse plane intersecting the center of the nuclear target 
     * (along the beamline direction, i.e. z). We project the horizontal
     * position from HP875 to the nuclear target since HPTG1 and HPTG2 are close
     * to the front of the target. We project the vertical position from VP875
     * since both VP873 and VP875 are sufficiently far from the nuclear target
     * and VP875 is closer to the nuclear target.
     * 
     * Order of horizontal BPMs: HP875, HPTG1, HPTG2
     * Order of vertical BPMs: VP873, VP875
     */
    Double_t horAng = (HPTGXOffset - HP875Offset) / (HPTGXZPos - HP875ZPos);
    horAng = atan(horAng);
    Double_t tgtHorPos = HP875Offset + horAng * (targetCenterZPos - HP875ZPos);

    Double_t verAng = (VP875Offset - VP873Offset) / (VP875ZPos - VP873ZPos);
    verAng = atan(verAng);
    Double_t tgtVerPos = VP873Offset + horAng * (targetCenterZPos - HP875ZPos);

    /**
     * Prioritize using MW876 data if it has reasonable fit parameters.
     * Otherwise, prioritize using MW875 data. If none of the multi-wire 
     * devices have reasonable fit parameters, default to using calcFoM()
     * instead of calcFoM2(). Note calcFoM() does not make use of any 
     * multi-wire data.
     * 
     * @note What counts as "reasonable" is up to interpretation. For the time
     * being, we use the hard-coded values that MicroBooNE used ( @see
     * getFOM2.cxx in MicroBooNE's ubraw GitHub repository for details).
     */

    std::optional<Double_t> tgtHorWidth = std::nullopt, tgtVerWidth = std::nullopt;
    if ( isReasonable( MW876HorWidth, 0.5, 10) 
        && isReasonable( MW876VerWidth, 0.3, 10) ) {
            tgtHorWidth = expandWidth( MW876HorWidth, p876x, 2);
            tgtVerWidth = expandWidth( MW876VerWidth, p876y, 2);

            /** @warning 
             * If specified N for expansion is greater than the number of
             * available coefficents, cannot calculate Figure of Merit
             */
            if ( !tgtHorWidth.has_value() || !tgtVerWidth.has_value()) {
                return -1.;
            }
    }
    else if ( isReasonable( MW875HorWidth, 0.5, 10) 
        && isReasonable(MW875VerWidth, 0.3, 10) ) {
            tgtHorWidth = expandWidth( MW875HorWidth, p875x, 2);
            tgtVerWidth = expandWidth( MW875VerWidth, p875y, 2);

            /** @warning 
             * If specified N for expansion is greater than the number of
             * available coefficents, cannot calculate Figure of Merit
             */
            if ( !tgtHorWidth.has_value() || !tgtVerWidth.has_value()) {
                return -1.;
            }
    }
    else { 
        //! @warning defaulting to Figure of Merit without multi-wire data
        std::cerr << "[WARNING] No reliable multi-wire data, defaulting to calcFoM()." << std::endl;
        return calcFoM( tgtHorPos, horAng, tgtVerPos, verAng, TOR);
    }

    //! Extract values from std::optional objects
    Double_t tgtHorWidthVal = tgtHorWidth.value();
    Double_t tgtVerWidthVal = tgtVerWidth.value();

    return calcFoM2( tgtHorPos, horAng, tgtVerPos, verAng, 
        TOR, tgtHorWidthVal, tgtVerWidthVal);
}









//!-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
//! @todo

/**
 * @brief Fit (hopefully) gaussian profile of a given multi-wire device and
 * extract the gaussian fit parameters.
 * 
 * @note This function should only be implemented if it is found that the
 * multi-wire fit parameters on IFBeam are not reliable. There's no sense in
 * re-fitting if we trust the IFBeam fit parameters. The IFBeam fit parameters
 * are what is referenced as kSpillMW875HorWidth, kSpillMW876VerPos, etc., in
 * the BNBVars.cxx and BNBVars.h files.
 */
void processMWRProfile() {
    std::cout << "Placeholder: extract gaussian fit params from MWR device \n";
}
