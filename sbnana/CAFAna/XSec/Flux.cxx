#include "sbnana/CAFAna/XSec/Flux.h"

#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Weight.h"

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
    // The issue with the CC interaction is there is a threshold around the
    // muon mass. This process works for lower neutrino energies too
    return nu->pdg == pdg && !nu->iscc &&
      nu->genie_mode == caf::kQE && nu->genie_inttype == caf::kNCQE &&
      nu->targetPDG == 1000180400 /* Argon 40 */ &&
      !nu->ischarm &&
      nu->hitnuc == 2112;
  }

  //----------------------------------------------------------------------
  const NuTruthWeight kInvXSec([](const caf::SRTrueInteractionProxy* nu)
                             {
                               // GENIE uses GeV internally. We ultimately want
                               // a flux in m^-2
                               const double GeV2perm2 = 2.56819e31;

                               return GeV2perm2/nu->xsec;
                             });

  //----------------------------------------------------------------------
  NuTruthCut IsNCQEOnArgonCut(int pdg)
  {
    return NuTruthCut([pdg](const caf::SRTrueInteractionProxy* nu)
                      {
                        return IsNCQEOnArgon(nu, pdg);
                      });
  }

  //----------------------------------------------------------------------
  FluxTimesNuclei::FluxTimesNuclei(INuTruthSource& src,
                                   const Binning& bins,
                                   const NuTruthCut& fidvol,
                                   int pdg,
                                   const NuTruthWeight& wgt)
    : Spectrum(src[IsNCQEOnArgonCut(pdg) && fidvol].Weighted(wgt *kInvXSec),
               NuTruthHistAxis("True neutrino energy (GeV)",
                               bins,
                               SIMPLENUTRUTHVAR(E))),
      fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  TH1D* FluxTimesNuclei::ToTH1(double pot,
                               Color_t col,
                               Style_t style,
                               EBinType bintype)
  {
    TH1D* ret = Spectrum::ToTH1(pot, col, style, kPOT, bintype);

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

  //----------------------------------------------------------------------
  EnsembleFluxTimesNuclei::EnsembleFluxTimesNuclei(INuTruthEnsembleSource& src,
                                                   const Binning& bins,
                                                   const NuTruthCut& fidvol,
                                                   int pdg,
                                                   const NuTruthWeight& wgt)
    : EnsembleSpectrum(src[IsNCQEOnArgonCut(pdg) && fidvol].Weighted(wgt * kInvXSec),
                       NuTruthHistAxis("True neutrino energy (GeV)",
                                       bins,
                                       SIMPLENUTRUTHVAR(E))),
      fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  TH1D* EnsembleFluxTimesNuclei::ToTH1(double pot,
                                       Color_t col,
                                       Style_t style,
                                       EBinType bintype)
  {
    TH1D* ret = EnsembleSpectrum::Nominal().ToTH1(pot, col, style, kPOT, bintype);

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
