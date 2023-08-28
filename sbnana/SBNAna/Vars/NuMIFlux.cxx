#include "sbnana/SBNAna/Vars/NuMIFlux.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ana {

  //// ----------------------------------------------
  NuMIPpfxFluxWeight::NuMIPpfxFluxWeight()
  {
    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIPpfxFluxWeight: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort();
    }

    fFluxFilePath = std::string(sbndata) +
                   "beamData/NuMIdata/2023-07-31_out_450.37_7991.98_79512.66_QEL11.root";

    TFile f(fFluxFilePath.c_str());
    if (f.IsZombie()) {
      std::cout << "NuMIPpfxFluxWeight: Failed to open " << fFluxFilePath << std::endl;
      std::abort();
    }

    for (int hcIdx : {0, 1}) {
      for (int flavIdx : {0, 1}) {
        for (int signIdx : {0, 1}) {
          std::string hNamePPFX = "ppfx_flux_weights/hweights_";
          if (hcIdx == 0)
            hNamePPFX += "fhc_";
          else
            hNamePPFX += "rhc_";
          if (flavIdx == 0)
            hNamePPFX += "nue";
          else
            hNamePPFX += "numu";
          if (signIdx == 1) hNamePPFX += "bar";

          TH1* h_ppfx = (TH1*)f.Get(hNamePPFX.c_str());
          if (!h_ppfx) {
            std::cout << "NuMIPpfxFluxWeight: failed to find " << hNamePPFX << " in " << f.GetName()
                      << std::endl;
            std::abort();
          }
          h_ppfx = (TH1*)h_ppfx->Clone(UniqueName().c_str());
          h_ppfx->SetDirectory(0);

          fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;
        }
      }
    }
  }

  double NuMIPpfxFluxWeight::GetNuWeight(const caf::Proxy<caf::SRTrueInteraction>& true_int) const {

    if ( abs(true_int.initpdg) == 16 ) return 1.0;

    if ( !fWeight[0][0][0] ) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = ( abs(true_int.initpdg) == 12 ) ? 0 : 1;
    unsigned int signIdx = ( true_int.initpdg > 0 ) ? 0 : 1;

    TH1* h = fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin( true_int.E );
    if(bin == 0 || bin == h->GetNbinsX()+1) return 1.0;
    return h->GetBinContent(bin);

  }

  NuMIPpfxFluxWeight::~NuMIPpfxFluxWeight()
  {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          delete fWeight[i][j][k];
  }

  //// ----------------------------------------------

  const Var kGetNuMIFluxWeight([](const caf::SRSliceProxy* slc) -> double {
    if (slc->truth.index < 0 || abs(slc->truth.initpdg) == 16) return 1.0;

    if (!FluxWeightNuMI.fWeight[0][0][0]) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = (abs(slc->truth.initpdg) == 12) ? 0 : 1;
    unsigned int signIdx = (slc->truth.initpdg > 0) ? 0 : 1;

    TH1* h = FluxWeightNuMI.fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(slc->truth.E);

    if (bin == 0 || bin == h->GetNbinsX() + 1) return 1.0;
    if (std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin))) return 1.0;
    return h->GetBinContent(bin);
  });

  //// ----------------------------------------------

}
