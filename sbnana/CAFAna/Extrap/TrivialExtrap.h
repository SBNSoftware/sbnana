#pragma once

#include "sbnana/CAFAna/Extrap/IExtrap.h"

#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/Weight.h"

namespace ana
{
  class Loaders;

  /// "Extrapolation" that simply returns the MC prediction
  class TrivialExtrap: public IExtrap
  {
  public:
    TrivialExtrap(ISliceSource& nonswapSrc,
                  ISliceSource& nueSrc,
                  ISliceSource& tauSrc,
                  ISliceSource& intrinsicSrc,
                  const HistAxis& axis);

    TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                  SpectrumLoaderBase& loaderNue,
                  SpectrumLoaderBase& loaderNuTau,
                  SpectrumLoaderBase& loaderIntrinsic,
                  const HistAxis& axis,
                  const Cut& cut,
                  const SystShifts& shift,
                  const Weight& wei);

    TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                  SpectrumLoaderBase& loaderNue,
                  SpectrumLoaderBase& loaderNuTau,
                  SpectrumLoaderBase& loaderIntrinsic,
                  const HistAxis& axis,
                  const SpillCut& spillcut,
                  const Cut& cut,
                  const SystShifts& shift,
                  const Weight& wei)
      : TrivialExtrap()
    {
      abort(); // TODO TODO TODO
    }

    TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                  SpectrumLoaderBase& loaderNue,
                  SpectrumLoaderBase& loaderNuTau,
                  SpectrumLoaderBase& loaderIntrinsic,
                  std::string label,
                  const Binning& bins,
                  const Var& var,
                  const Cut& cut,
                  const SystShifts& shift,
                  const Weight& wei);

    TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                  SpectrumLoaderBase& loaderNue,
                  SpectrumLoaderBase& loaderNuTau,
                  SpectrumLoaderBase& loaderIntrinsic,
                  std::string label,
                  const Binning& bins,
                  const Var& var,
                  const SpillCut& spillcut,
                  const Cut& cut,
                  const SystShifts& shift,
                  const Weight& wei)
      : TrivialExtrap()
    {
      abort(); // TODO TODO TODO
    }

    TrivialExtrap(Loaders& loaders,
                  std::string label,
                  const Binning& bins,
                  const Var& var,
                  const Cut& cut,
                  const SystShifts& shift = kNoShift,
                  const Weight& wei = kUnweighted);

    TrivialExtrap(Loaders& loaders,
                  std::string label,
                  const Binning& bins,
                  const Var& var,
                  const SpillCut& spillcut,
                  const Cut& cut,
                  const SystShifts& shift = kNoShift,
                  const Weight& wei = kUnweighted)
      : TrivialExtrap()
    {
      abort(); // TODO TODO TODO
    }

    TrivialExtrap(Loaders& loaders,
                  const HistAxis& axis,
                  const Cut& cut,
                  const SystShifts& shift = kNoShift,
                  const Weight& wei = kUnweighted);

    TrivialExtrap(Loaders& loaders,
                  const HistAxis& axis,
                  const SpillCut& spillcut,
                  const Cut& cut,
                  const SystShifts& shift = kNoShift,
                  const Weight& wei = kUnweighted)
      : TrivialExtrap()
    {
      abort(); // TODO TODO TODO
    }

    virtual OscillatableSpectrum NueSurvComponent()       {return fNueSurv;}
    virtual OscillatableSpectrum AntiNueSurvComponent()   {return fNueSurvAnti;}

    virtual OscillatableSpectrum NumuSurvComponent()      {return fNumuSurv;}
    virtual OscillatableSpectrum AntiNumuSurvComponent()  {return fNumuSurvAnti;}

    virtual OscillatableSpectrum NueAppComponent()        {return fNueApp;}
    virtual OscillatableSpectrum AntiNueAppComponent()    {return fNueAppAnti;}

    virtual OscillatableSpectrum NumuAppComponent()       {return fNumuApp;}
    virtual OscillatableSpectrum AntiNumuAppComponent()   {return fNumuAppAnti;}

    virtual OscillatableSpectrum TauFromEComponent()      {return fTauFromE;}
    virtual OscillatableSpectrum AntiTauFromEComponent()  {return fTauFromEAnti;}

    virtual OscillatableSpectrum TauFromMuComponent()     {return fTauFromMu;}
    virtual OscillatableSpectrum AntiTauFromMuComponent() {return fTauFromMuAnti;}

    virtual OscillatableSpectrum NCComponentFromNumu() {return fNCFromNumu;}
    virtual OscillatableSpectrum NCComponentFromNue() {return fNCFromNue;}

    virtual void SaveTo(TDirectory* dir) const;
    static std::unique_ptr<TrivialExtrap> LoadFrom(TDirectory* dir);

  protected:
    TrivialExtrap()
      : fNueApp(OscillatableSpectrum::Uninitialized()),    fNueAppAnti(OscillatableSpectrum::Uninitialized()),
        fNumuSurv(OscillatableSpectrum::Uninitialized()),  fNumuSurvAnti(OscillatableSpectrum::Uninitialized()),
        fNumuApp(OscillatableSpectrum::Uninitialized()),   fNumuAppAnti(OscillatableSpectrum::Uninitialized()),
        fNueSurv(OscillatableSpectrum::Uninitialized()),   fNueSurvAnti(OscillatableSpectrum::Uninitialized()),
        fTauFromE(OscillatableSpectrum::Uninitialized()),  fTauFromEAnti(OscillatableSpectrum::Uninitialized()),
        fTauFromMu(OscillatableSpectrum::Uninitialized()), fTauFromMuAnti(OscillatableSpectrum::Uninitialized()),
        fNCFromNumu(OscillatableSpectrum::Uninitialized()), fNCFromNue(OscillatableSpectrum::Uninitialized())
    {
    }

    OscillatableSpectrum fNueApp,    fNueAppAnti;
    OscillatableSpectrum fNumuSurv,  fNumuSurvAnti;
    OscillatableSpectrum fNumuApp,   fNumuAppAnti;
    OscillatableSpectrum fNueSurv,   fNueSurvAnti;
    OscillatableSpectrum fTauFromE,  fTauFromEAnti;
    OscillatableSpectrum fTauFromMu, fTauFromMuAnti;
    OscillatableSpectrum fNCFromNumu, fNCFromNue;
  };
}
