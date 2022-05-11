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
  PredictionNoOsc::PredictionNoOsc(ISliceSource& src, const HistAxis& axis)
    : fSpectrum       (src, axis),
      fSpectrumNC     (src[ kIsNC], axis),
      fSpectrumNumu   (src[!kIsNC][kIsNumuCC ][!kIsAntiNu], axis),
      fSpectrumNumubar(src[!kIsNC][kIsNumuCC ][ kIsAntiNu], axis),
      fSpectrumNue    (src[!kIsNC][kIsBeamNue][!kIsAntiNu], axis),
      fSpectrumNuebar (src[!kIsNC][kIsBeamNue][ kIsAntiNu], axis)
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
  void PredictionNoOsc::SaveTo(TDirectory* dir, const std::string& name) const
  {
    TDirectory* tmp = gDirectory;

    dir = dir->mkdir(name.c_str()); // switch to subdir
    dir->cd();

    TObjString("PredictionNoOsc").Write("type");

    fSpectrum.SaveTo(dir, "spect");
    fSpectrumNC.SaveTo(dir, "spect_nc");
    fSpectrumNumu.SaveTo(dir, "spect_numu");
    fSpectrumNumubar.SaveTo(dir, "spect_numubar");
    fSpectrumNue.SaveTo(dir, "spect_nue");
    fSpectrumNuebar.SaveTo(dir, "spect_nuebar");

    dir->Write();
    delete dir;

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<PredictionNoOsc> PredictionNoOsc::LoadFrom(TDirectory* dir, const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str()); // switch to subdir
    assert(dir);

    PredictionNoOsc* ret = new PredictionNoOsc(
      *ana::LoadFrom<Spectrum>(dir, "spect"),
      *ana::LoadFrom<Spectrum>(dir, "spect_nc"),
      *ana::LoadFrom<Spectrum>(dir, "spect_numu"),
      *ana::LoadFrom<Spectrum>(dir, "spect_numubar"),
      *ana::LoadFrom<Spectrum>(dir, "spect_nue"),
      *ana::LoadFrom<Spectrum>(dir, "spect_nuebar"));

    // Can't use make_unique because constructor is protected
    return std::unique_ptr<PredictionNoOsc>(ret);
  }
}
