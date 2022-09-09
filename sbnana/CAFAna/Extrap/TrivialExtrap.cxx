#include "sbnana/CAFAna/Extrap/TrivialExtrap.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  // TODO it's not ideal this has to be duplicated for a third time
  namespace{
    const Var kTrueE([](const caf::SRSliceProxy* slc) -> double {return slc->truth.E;});

    const Var kBaseline([](const caf::SRSliceProxy* slc) -> double {return slc->truth.baseline * 1e-3;}); // m -> km

    const Var kTrueLOverE = kBaseline / kTrueE;
  }

  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(ISliceSource& nonswapSrc,
                               ISliceSource& nueSrc,
                               ISliceSource& tauSrc,
                               ISliceSource& intrinsicSrc,
                               const HistAxis& axis)
    :
    fNueApp       (nueSrc      [kIsNueApp    && !kIsAntiNu], axis),
    fNueAppAnti   (nueSrc      [kIsNueApp    &&  kIsAntiNu], axis),
    fNumuSurv     (nonswapSrc  [kIsNumuCC    && !kIsAntiNu], axis),
    fNumuSurvAnti (nonswapSrc  [kIsNumuCC    &&  kIsAntiNu], axis),
    fNumuApp      (tauSrc      [kIsNumuApp   && !kIsAntiNu], axis),
    fNumuAppAnti  (tauSrc      [kIsNumuApp   &&  kIsAntiNu], axis),
    fNueSurv      (intrinsicSrc[kIsBeamNue   && !kIsAntiNu], axis),
    fNueSurvAnti  (intrinsicSrc[kIsBeamNue   &&  kIsAntiNu], axis),
    fTauFromE     (nueSrc      [kIsTauFromE  && !kIsAntiNu], axis),
    fTauFromEAnti (nueSrc      [kIsTauFromE  &&  kIsAntiNu], axis),
    fTauFromMu    (tauSrc      [kIsTauFromMu && !kIsAntiNu], axis),
    fTauFromMuAnti(tauSrc      [kIsTauFromMu &&  kIsAntiNu], axis),
    fNCFromNumu   (nonswapSrc[kIsNCFromNumu], axis),
    fNCFromNue    (nonswapSrc[kIsNCFromNue ], axis)
  {
    // All swapped files are equally valid as a source of NCs. This
    // approximately doubles/triples our statistics. SpectrumLoader just adds
    // events and POT for both cases, which is the right thing to do.

    nueSrc[kIsNCFromNumu].GetVars(axis.GetVar1D(), kTrueLOverE).Register(&fNCFromNumu);
    tauSrc[kIsNCFromNumu].GetVars(axis.GetVar1D(), kTrueLOverE).Register(&fNCFromNumu);

    nueSrc[kIsNCFromNue].GetVars(axis.GetVar1D(), kTrueLOverE).Register(&fNCFromNue);
    tauSrc[kIsNCFromNue].GetVars(axis.GetVar1D(), kTrueLOverE).Register(&fNCFromNue);

    //Also load in intrinsic nues from nonswap file
    nonswapSrc[kIsBeamNue && !kIsAntiNu].GetVars(axis.GetVar1D(), kTrueLOverE).Register(&fNueSurv);
    nonswapSrc[kIsBeamNue &&  kIsAntiNu].GetVars(axis.GetVar1D(), kTrueLOverE).Register(&fNueSurvAnti);
  }

  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(SliceSources& srcs, const HistAxis& axis)
    : TrivialExtrap(srcs.GetLoader(kMC, kNonSwap),
                    srcs.GetLoader(kMC, kNueSwap),
                    srcs.GetLoader(kMC, kNuTauSwap),
		    srcs.GetLoader(kMC, kIntrinsic),
                    axis)
  {
  }

  //----------------------------------------------------------------------
  void TrivialExtrap::SaveTo(TDirectory* dir, const std::string& name) const
  {
    TDirectory* tmp = gDirectory;

    dir = dir->mkdir(name.c_str()); // switch to subdir
    dir->cd();

    TObjString("TrivialExtrap").Write("type");

    fNueApp.SaveTo(dir, "nue_app");
    fNueAppAnti.SaveTo(dir, "nue_app_anti");
    fNCFromNumu.SaveTo(dir, "nc_from_numu");
    fNCFromNue.SaveTo(dir, "nc_from_nue");
    fNumuSurv.SaveTo(dir, "numu_surv");
    fNumuSurvAnti.SaveTo(dir, "numu_surv_anti");
    fNumuApp.SaveTo(dir, "numu_app");
    fNumuAppAnti.SaveTo(dir, "numu_app_anti");
    fNueSurv.SaveTo(dir, "nue_surv");
    fNueSurvAnti.SaveTo(dir, "nue_surv_anti");
    fTauFromE.SaveTo(dir, "nutau_from_nue");
    fTauFromEAnti.SaveTo(dir, "nutau_from_nue_anti");
    fTauFromMu.SaveTo(dir, "nutau_from_numu");
    fTauFromMuAnti.SaveTo(dir, "nutau_from_numu_anti");

    dir->Write();
    delete dir;

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<TrivialExtrap> TrivialExtrap::LoadFrom(TDirectory* dir, const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str()); // switch to subdir
    assert(dir);

    std::unique_ptr<TrivialExtrap> ret(new TrivialExtrap);

    ret->fNueApp        = *OscillatableSpectrum::LoadFrom(dir, "nue_app");
    ret->fNueAppAnti    = *OscillatableSpectrum::LoadFrom(dir, "nue_app_anti");
    ret->fNumuSurv      = *OscillatableSpectrum::LoadFrom(dir, "numu_surv");
    ret->fNumuSurvAnti  = *OscillatableSpectrum::LoadFrom(dir, "numu_surv_anti");
    ret->fNumuApp       = *OscillatableSpectrum::LoadFrom(dir, "numu_app");
    ret->fNumuAppAnti   = *OscillatableSpectrum::LoadFrom(dir, "numu_app_anti");
    ret->fNueSurv       = *OscillatableSpectrum::LoadFrom(dir, "nue_surv");
    ret->fNueSurvAnti   = *OscillatableSpectrum::LoadFrom(dir, "nue_surv_anti");
    ret->fTauFromE      = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_nue");
    ret->fTauFromEAnti  = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_nue_anti");
    ret->fTauFromMu     = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_numu");
    ret->fTauFromMuAnti = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_numu_anti");
    ret->fNCFromNumu    = *OscillatableSpectrum::LoadFrom(dir, "nc_from_numu");
    ret->fNCFromNue     = *OscillatableSpectrum::LoadFrom(dir, "nc_from_nue");

    return ret;
  }
}
