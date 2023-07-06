// BH - 2023
// HEAVILY based on NuMI flux syst, and thanks to Tony Wood for discussing the right histograms to use

#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include <string_view>


class TH1;

namespace ana
{

  class NuMIPpfxFluxWeight
  {
    static constexpr std::string_view fluxFileName = "2023-07-06_out_450.37_7991.98_79512.66_QEL11.root";

  public:
    NuMIPpfxFluxWeight();
    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
  };

  // set up to use the flux weight
  NuMIPpfxFluxWeight FluxWeightNuMI;
  extern const Var kGetNuMIFluxWeight;

}
