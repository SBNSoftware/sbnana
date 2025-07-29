// BH - 2023
// HEAVILY based on NuMI flux syst, and thanks to Tony Wood for discussing the right histograms to use

#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include <string>

class TH1;

namespace ana
{

  //===============================================================
  // 05/19/25 JK) This one uses "2025-04-08_out_450.37_7991.98_79512.66.root",
  //              to provide a CV correction of beam width setting 1.5mm/1.4mm;
  //              NuMI2023 reprocessing simulation was done with CV beam width of 1.4mm,
  //              but the data is more like 1.5mm.
  class NuMIFluxCorrection
  {
  public:
    NuMIFluxCorrection();
    ~NuMIFluxCorrection();

    double GetWeightFromSRTrueInt(const caf::SRTrueInteractionProxy* nu) const;
    unsigned int ParentPDGToIdx(int pdg) const;

    // PPFX correction
    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
    // Additional CV correction from NuMI reproc-to-PPFXCalculationNominal of this file
    mutable TH1* fWeightCVCorr[2][2][2][4]; // [fhc/rhc][nue/numu][nu/nubar][parent pid (pipm/kpm/k0l/mu)]

    static NuMIFluxCorrection& Instance();

  protected:
    std::string fFluxFilePath;
  };

  extern const Var kGetNuMIFluxCorrection;
  extern const TruthVar kGetTruthNuMIFluxCorrection;


}
