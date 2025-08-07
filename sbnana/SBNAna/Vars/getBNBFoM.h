//! ////////////////////////////////////////////////////////////////////////////
//! @file: getBNBFoM.h                                                           
//! @author: Jacob Smith (smithja)
//! @email: jacob.a.smith@stonybrook.edu                                       
//!                                                                           
//! @brief Header file to define all of functions used in calculating a Booster
//! Neutrino Beam (BNB) Figure of Merit (FoM) for use on the SBN experiments.                  
//! ////////////////////////////////////////////////////////////////////////////
#ifndef _GETBNBFOM_H
#define _GETBNBFOM_H

#include <optional>

#include "TROOT.h"

std::optional<Double_t> getValidBPMCalibVal( 
    const int& timestamp, 
    const std::string& device
);

std::optional<Double_t> getValidBPMCalibTime( 
    const int& timestamp, 
    const std::string& device
);
    
void swimBNB(
    const Double_t centroid1[6], 
    const Double_t sigma1[6][6],
    const Double_t xferc[6][6], 
    const Double_t xfers[6][6],
    Double_t& cx,
    Double_t& cy, 
    Double_t& sx, 
    Double_t& sy, 
    Double_t& rho
);

Double_t funcIntbivar(
    const Double_t cx, 
    const Double_t cy,
    const Double_t sx, 
    const Double_t sy,
    const Double_t rho
);

Double_t calcFoM(
    const Double_t horPos, 
    const Double_t horAng, 
    const Double_t verPos, 
    const Double_t verAng, 
    const Double_t PPP
);

Double_t calcFoM2(
    const Double_t horPos, 
    const Double_t horAng, 
    const Double_t verPos, 
    const Double_t verAng, 
    const Double_t PPP, 
    const Double_t tgtHorSigma, 
    const Double_t tgtVerSigma
);

Double_t getBNBFoM( 
    const Double_t spillTimeSec,
    const Double_t TOR860,
    const Double_t TOR875,
    const Double_t HP875, 
    const Double_t HPTG1,
    const Double_t HPTG2,
    const Double_t VP873,
    const Double_t VP875
);

Double_t getBNBFoM2( 
    const Double_t spillTimeSec,
    const Double_t TOR860,
    const Double_t TOR875,
    const Double_t HP875,
    const Double_t HPTG1, 
    const Double_t HPTG2,
    const Double_t VP873,
    const Double_t VP875,
    const Double_t MW875HorSigma,
    const Double_t MW875VerSigma,
    const Double_t MW876HorSigma,
    const Double_t MW876VerSigma
);

#endif