#include "sbnana/CAFAna/XSec/Flux.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TH1.h"

namespace ana
{
  //----------------------------------------------------------------------
  bool IsCCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg)
  {
    // It doesn't matter exactly what process we choose, so long as it
    // corresponds to one of the modes GENIE considers as a top-level
    // cross-section. Here we go for CC QE on Argon.
    return nu->pdg == pdg && nu->iscc &&
      nu->genie_mode == caf::kQE && nu->genie_inttype == caf::kCCQE &&
      nu->targetPDG == 1000180400 /* Argon 40 */ &&
      !nu->ischarm;
  }

  //----------------------------------------------------------------------
  bool IsNCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg)
  {
    // This process works for very low neutrino energy where the muon mass
    // becomes relevant.
    return nu->pdg == pdg && !nu->iscc &&
      nu->genie_mode == caf::kQE && nu->genie_inttype == caf::kNCQE &&
      nu->targetPDG == 1000180400 /* Argon 40 */ && 
      !nu->ischarm &&
      nu->hitnuc == 2112;
  }

  //----------------------------------------------------------------------
  const Var kInvXSec([](const caf::SRSliceProxy* sr)
                     {
                       // GENIE uses GeV internally. We ultimately want a flux
                       // in m^-2
                       const double GeV2perm2 = 2.56819e31;

                       return GeV2perm2/sr->truth.xsec;
                     });

  //----------------------------------------------------------------------
  // TODO can this operate completely in true interaction mode? Right now we
  // are folding in a slicing efficiency
  Cut IsCCQEOnArgonCut(int pdg)
  {
    return Cut([pdg](const caf::SRSliceProxy* slc)
               {
                 if(slc->truth.index < 0) return false;
                 return IsCCQEOnArgon(&slc->truth, pdg);
               });
  }

  //----------------------------------------------------------------------
  FluxTimesNuclei::FluxTimesNuclei(SpectrumLoaderBase& loader,
                                   const Binning& bins,
                                   const Cut& fidvol,
                                   int pdg)
    : Spectrum("", bins, loader, SIMPLEVAR(truth.E),
               SIMPLEVAR(truth.index) >= 0 && IsCCQEOnArgonCut(pdg) && fidvol,
               kNoShift, kInvXSec),
      fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  TH1D* FluxTimesNuclei::ToTH1(double exposure,
                               Color_t col,
                               Style_t style,
                               EBinType bintype)
  {
    TH1D* ret = Spectrum::ToTH1(exposure, col, style, kPOT, bintype);
    ret->GetXaxis()->SetTitle("True neutrino energy (GeV)");

    std::string ytitle = "Flux #times nuclei (";
    switch(fPdg){
    case +12: ytitle += "#nu_e"; break;
    case -12: ytitle += "#bar{#nu}_e"; break;
    case +14: ytitle += "#nu_{#mu}"; break;
    case -14: ytitle += "#bar{#nu}_{#mu}"; break;
    }

    ytitle += "/ m^{2}";
    if(bintype == kBinDensity) ytitle += " / GeV";
    ytitle += ")";
    ret->GetYaxis()->SetTitle(ytitle.c_str());
    return ret;
  }
}
