#include "sbnana/CAFAna/Extrap/TrivialExtrap.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                               SpectrumLoaderBase& loaderNue,
                               SpectrumLoaderBase& loaderNuTau,
                               SpectrumLoaderBase& loaderIntrinsic,
                               const HistAxis& axis,
                               const SpillCut& spillcut,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    :
    fNueApp       (loaderNue,     axis, spillcut, cut && kIsNueApp    && !kIsAntiNu, shift, wei),
    fNueAppAnti   (loaderNue,     axis, spillcut, cut && kIsNueApp    &&  kIsAntiNu, shift, wei),

    fNumuSurv     (loaderNonswap, axis, spillcut, cut && kIsNumuCC    && !kIsAntiNu, shift, wei),
    fNumuSurvAnti (loaderNonswap, axis, spillcut, cut && kIsNumuCC    &&  kIsAntiNu, shift, wei),

    fNumuApp      (loaderNuTau,   axis, spillcut, cut && kIsNumuApp   && !kIsAntiNu, shift, wei),
    fNumuAppAnti  (loaderNuTau,   axis, spillcut, cut && kIsNumuApp   &&  kIsAntiNu, shift, wei),

    fNueSurv      (loaderIntrinsic, axis, spillcut, cut && kIsBeamNue   && !kIsAntiNu, shift, wei),
    fNueSurvAnti  (loaderIntrinsic, axis, spillcut, cut && kIsBeamNue   &&  kIsAntiNu, shift, wei),

    fTauFromE     (loaderNue,     axis, spillcut, cut && kIsTauFromE  && !kIsAntiNu, shift, wei),
    fTauFromEAnti (loaderNue,     axis, spillcut, cut && kIsTauFromE  &&  kIsAntiNu, shift, wei),

    fTauFromMu    (loaderNuTau,   axis, spillcut, cut && kIsTauFromMu && !kIsAntiNu, shift, wei),
    fTauFromMuAnti(loaderNuTau,   axis, spillcut, cut && kIsTauFromMu &&  kIsAntiNu, shift, wei),

    fNCFromNumu   (loaderNonswap, axis, spillcut, cut && kIsNCFromNumu,     shift, wei),
    fNCFromNue    (loaderNonswap, axis, spillcut, cut && kIsNCFromNue,      shift, wei)
  {
    // All swapped files are equally valid as a source of NCs. This
    // approximately doubles/triples our statistics. SpectrumLoader just adds
    // events and POT for both cases, which is the right thing to do.

    loaderNue  .AddReweightableSpectrum(fNCFromNumu, axis.GetMultiDVar(), spillcut, cut && kIsNCFromNumu, shift, wei);
    loaderNuTau.AddReweightableSpectrum(fNCFromNumu, axis.GetMultiDVar(), spillcut, cut && kIsNCFromNumu, shift, wei);

    loaderNue  .AddReweightableSpectrum(fNCFromNue, axis.GetMultiDVar(), spillcut, cut && kIsNCFromNue, shift, wei);
    loaderNuTau.AddReweightableSpectrum(fNCFromNue, axis.GetMultiDVar(), spillcut, cut && kIsNCFromNue, shift, wei);

    //Also load in intrinsic nues from nonswap file
    loaderNonswap.AddReweightableSpectrum(fNueSurv, axis.GetMultiDVar(), spillcut, cut && kIsBeamNue && !kIsAntiNu, shift, wei);
    loaderNonswap.AddReweightableSpectrum(fNueSurvAnti, axis.GetMultiDVar(), spillcut, cut && kIsBeamNue && kIsAntiNu, shift, wei);

  }


  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                               SpectrumLoaderBase& loaderNue,
                               SpectrumLoaderBase& loaderNuTau,
                               SpectrumLoaderBase& loaderIntrinsic,
                               std::string label,
                               const Binning& bins,
                               const Var& var,
                               const SpillCut& spillcut,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    :
    TrivialExtrap(loaderNonswap, loaderNue, loaderNuTau, loaderIntrinsic,
                  HistAxis(label, bins, var),
                  spillcut, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(Loaders& loaders,
                               std::string label,
                               const Binning& bins,
                               const Var& var,
                               const SpillCut& spillcut,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    : TrivialExtrap(loaders, HistAxis(label, bins, var), spillcut, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(Loaders& loaders,
                               const HistAxis& axis,
                               const SpillCut& spillcut,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    : TrivialExtrap(loaders.GetLoader(Loaders::kMC, ana::kBeam, Loaders::kNonSwap),
                    loaders.GetLoader(Loaders::kMC, ana::kBeam, Loaders::kNueSwap),
                    loaders.GetLoader(Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap),
		    loaders.GetLoader(Loaders::kMC, ana::kBeam, Loaders::kIntrinsic),
                    axis, spillcut, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  void TrivialExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("TrivialExtrap").Write("type");

    fNueApp.SaveTo(dir->mkdir("nue_app"));
    fNueAppAnti.SaveTo(dir->mkdir("nue_app_anti"));
    fNCFromNumu.SaveTo(dir->mkdir("nc_from_numu"));
    fNCFromNue.SaveTo(dir->mkdir("nc_from_nue"));
    fNumuSurv.SaveTo(dir->mkdir("numu_surv"));
    fNumuSurvAnti.SaveTo(dir->mkdir("numu_surv_anti"));
    fNumuApp.SaveTo(dir->mkdir("numu_app"));
    fNumuAppAnti.SaveTo(dir->mkdir("numu_app_anti"));
    fNueSurv.SaveTo(dir->mkdir("nue_surv"));
    fNueSurvAnti.SaveTo(dir->mkdir("nue_surv_anti"));
    fTauFromE.SaveTo(dir->mkdir("nutau_from_nue"));
    fTauFromEAnti.SaveTo(dir->mkdir("nutau_from_nue_anti"));
    fTauFromMu.SaveTo(dir->mkdir("nutau_from_numu"));
    fTauFromMuAnti.SaveTo(dir->mkdir("nutau_from_numu_anti"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<TrivialExtrap> TrivialExtrap::LoadFrom(TDirectory* dir)
  {
    std::unique_ptr<TrivialExtrap> ret(new TrivialExtrap);

    // This is a lot of repetitive typing. Define some macros
#define LOAD_OSC(FIELD, LABEL) assert(dir->GetDirectory(LABEL)); ret->FIELD = *OscillatableSpectrum::LoadFrom(dir->GetDirectory(LABEL));
#define LOAD_SPECT(FIELD, LABEL) assert(dir->GetDirectory(LABEL)); ret->FIELD = *Spectrum::LoadFrom(dir->GetDirectory(LABEL));

    LOAD_OSC(fNueApp,        "nue_app");
    LOAD_OSC(fNueAppAnti,    "nue_app_anti");
    LOAD_OSC(fNumuSurv,      "numu_surv");
    LOAD_OSC(fNumuSurvAnti,  "numu_surv_anti");
    LOAD_OSC(fNumuApp,       "numu_app");
    LOAD_OSC(fNumuAppAnti,   "numu_app_anti");
    LOAD_OSC(fNueSurv,       "nue_surv");
    LOAD_OSC(fNueSurvAnti,   "nue_surv_anti");
    LOAD_OSC(fTauFromE,      "nutau_from_nue");
    LOAD_OSC(fTauFromEAnti,  "nutau_from_nue_anti");
    LOAD_OSC(fTauFromMu,     "nutau_from_numu");
    LOAD_OSC(fTauFromMuAnti, "nutau_from_numu_anti");
    LOAD_OSC(fNCFromNumu,    "nc_from_numu");
    LOAD_OSC(fNCFromNue,     "nc_from_nue");

    return ret;
  }
}
