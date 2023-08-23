// BH - 2023
// HEAVILY based on NuMI flux syst, and thanks to Tony Wood for discussing the right histograms to use

#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include <string>

class TH1;

namespace ana
{

  class NuMIPpfxFluxWeight
  {
  public:
    NuMIPpfxFluxWeight();
    ~NuMIPpfxFluxWeight();
    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
    double GetNuWeight(const caf::Proxy<caf::SRTrueInteraction>& true_int) const;

  protected:
    std::string fFluxFilePath;
  };

  // set up to use the flux weight
  static const NuMIPpfxFluxWeight FluxWeightNuMI;
  extern const Var kGetNuMIFluxWeight;

}
