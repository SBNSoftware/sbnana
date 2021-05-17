#include "sbnana/CAFAna/Prediction/PredictionNoOsc.h"

#include "sbnana/CAFAna/Extrap/IExtrap.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"

#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Extrap/TrivialExtrap.h"

#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "TDirectory.h"
#include "TObjString.h"


namespace ana
{
  //----------------------------------------------------------------------
  PredictionNoOsc::PredictionNoOsc(SpectrumLoaderBase& loader,
                                   const std::string& label,
                                   const Binning& bins,
                                   const Var& var,
                                   const SpillCut& spillcut,
                                   const Cut& cut,
                                   const SystShifts& shift,
                                   const Var& wei)
    : PredictionNoOsc(loader, HistAxis(label, bins, var), spillcut, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  PredictionNoOsc::PredictionNoOsc(SpectrumLoaderBase& loader,
                                   const HistAxis& axis,
                                   const SpillCut& spillcut,
                                   const Cut& cut,
                                   const SystShifts& shift,
                                   const Var& wei)
    : fSpectrum(       loader, axis, spillcut, cut,                                        shift, wei),
      fSpectrumNC(     loader, axis, spillcut, cut &&  kIsNC,                              shift, wei),
      fSpectrumNumu(   loader, axis, spillcut, cut && !kIsNC && kIsNumuCC &&  !kIsAntiNu,  shift, wei),
      fSpectrumNumubar(loader, axis, spillcut, cut && !kIsNC && kIsNumuCC &&   kIsAntiNu,  shift, wei),
      fSpectrumNue(    loader, axis, spillcut, cut && !kIsNC && kIsBeamNue && !kIsAntiNu,  shift, wei),
      fSpectrumNuebar( loader, axis, spillcut, cut && !kIsNC && kIsBeamNue &&  kIsAntiNu,  shift, wei)
  {
  }

  //----------------------------------------------------------------------
  Spectrum PredictionNoOsc::PredictComponent(osc::IOscCalc* /*calc*/,
                                             Flavors::Flavors_t flav,
                                             Current::Current_t curr,
                                             Sign::Sign_t sign) const
  {
    if(flav == Flavors::kAll &&
       curr == Current::kBoth &&
       sign == Sign::kBoth)
      return Predict(0); // Faster

    if(curr & Current::kNC){
      // We don't have NC broken down by sign or flavour
      assert(flav & Flavors::kAll && sign & Sign::kBoth);
      return fSpectrumNC;
    }

    assert(curr == Current::kCC);

    using namespace Flavors;
    using namespace Current;
    using namespace Sign;

    Spectrum ret = fSpectrum;
    ret.Clear();

    // Safe to assume by this point that it's 100% CC
    if(flav & kNuMuToNuMu && sign & kNu)     ret += fSpectrumNumu;
    if(flav & kNuMuToNuMu && sign & kAntiNu) ret += fSpectrumNumubar;
    if(flav & kNuEToNuE   && sign & kNu)     ret += fSpectrumNue;
    if(flav & kNuEToNuE   && sign & kAntiNu) ret += fSpectrumNuebar;

    return ret;
  }

  //----------------------------------------------------------------------
  void PredictionNoOsc::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionNoOsc").Write("type");

    fSpectrum.SaveTo(dir->mkdir("spect"));
    fSpectrumNC.SaveTo(dir->mkdir("spect_nc"));
    fSpectrumNumu.SaveTo(dir->mkdir("spect_numu"));
    fSpectrumNumubar.SaveTo(dir->mkdir("spect_numubar"));
    fSpectrumNue.SaveTo(dir->mkdir("spect_nue"));
    fSpectrumNuebar.SaveTo(dir->mkdir("spect_nuebar"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<PredictionNoOsc> PredictionNoOsc::LoadFrom(TDirectory* dir)
  {    
    PredictionNoOsc* ret = new PredictionNoOsc(
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("spect")),
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("spect_nc")),
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("spect_numu")),
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("spect_numubar")),
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("spect_nue")),
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("spect_nuebar")));

    // Can't use make_unique because constructor is protected
    return std::unique_ptr<PredictionNoOsc>(ret);
  }
}
