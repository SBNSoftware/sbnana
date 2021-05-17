#include "sbnana/CAFAna/Prediction/PredictionSBNExtrap.h"

#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Ratio.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1D.h"

namespace ana
{
  //----------------------------------------------------------------------
  PredictionSBNExtrap::PredictionSBNExtrap(Loaders& loadersND,
                                           Loaders& loadersFD,
                                           const HistAxis& axis,
                                           const SpillCut& spillcut,
                                           const Cut& cut,
                                           const SystShifts& shift_mc,
                                           const Var& wei_mc,
                                           const SystShifts& shift_data,
                                           const Var& wei_data)
    : fPredND(loadersND, axis, spillcut, cut, shift_mc, wei_mc),
      fPredFD(loadersFD, axis, spillcut, cut, shift_mc, wei_mc),
      // TODO how on earth do fake ND oscillations get in here?
      fDataND(loadersND.GetLoader(Loaders::kData), axis, spillcut, cut, shift_data, wei_data)
  {
  }

  //----------------------------------------------------------------------
  PredictionSBNExtrap::~PredictionSBNExtrap()
  {
  }

  //----------------------------------------------------------------------
  Spectrum PredictionSBNExtrap::Predict(osc::IOscCalc* calc) const
  {
    return PredictComponent(calc,
                            Flavors::kAll,
                            Current::kBoth,
                            Sign::kBoth);
  }

  //----------------------------------------------------------------------
  Spectrum PredictionSBNExtrap::PredictComponent(osc::IOscCalc* calc,
                                                 Flavors::Flavors_t flav,
                                                 Current::Current_t curr,
                                                 Sign::Sign_t sign) const
  {
    ((osc::IOscCalcAdjustable*)calc)->SetL(kBaselineSBND);
    const Ratio r = fDataND / fPredND.Predict(calc);
    ((osc::IOscCalcAdjustable*)calc)->SetL(kBaselineIcarus);
    return r * fPredFD.PredictComponent(calc, flav, curr, sign);
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum PredictionSBNExtrap::ComponentCC(int from, int to) const
  {
    /*
    if(from == +12 && to == +12) return fSBNExtrap->NueSurvComponent();
    if(from == -12 && to == -12) return fSBNExtrap->AntiNueSurvComponent();

    if(from == +12 && to == +14) return fSBNExtrap->NumuAppComponent();
    if(from == -12 && to == -14) return fSBNExtrap->AntiNumuAppComponent();

    if(from == +12 && to == +16) return fSBNExtrap->TauFromEComponent();
    if(from == -12 && to == -16) return fSBNExtrap->AntiTauFromEComponent();

    if(from == +14 && to == +12) return fSBNExtrap->NueAppComponent();
    if(from == -14 && to == -12) return fSBNExtrap->AntiNueAppComponent();

    if(from == +14 && to == +14) return fSBNExtrap->NumuSurvComponent();
    if(from == -14 && to == -14) return fSBNExtrap->AntiNumuSurvComponent();

    if(from == +14 && to == +16) return fSBNExtrap->TauFromMuComponent();
    if(from == -14 && to == -16) return fSBNExtrap->AntiTauFromMuComponent();
    */
    assert(0 && "Not reached");
  }

  //----------------------------------------------------------------------
  // Spectrum PredictionSBNExtrap::ComponentNC() const
  // {
  //   return fSBNExtrap->NCComponent();
  // }

  //----------------------------------------------------------------------
  void PredictionSBNExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionSBNExtrap").Write("type");

    fPredND.SaveTo(dir->mkdir("predND"));
    fPredFD.SaveTo(dir->mkdir("predFD"));
    fDataND.SaveTo(dir->mkdir("dataND"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<PredictionSBNExtrap> PredictionSBNExtrap::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("predND"));
    assert(dir->GetDirectory("predFD"));
    assert(dir->GetDirectory("dataND"));

    // TODO are these leaks?
    return std::unique_ptr<PredictionSBNExtrap>(new PredictionSBNExtrap(
      *ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("predND")).release(),
      *ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("predFD")).release(),
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("dataND")).release()));
  }

  //----------------------------------------------------------------------
  SBNExtrapGenerator::SBNExtrapGenerator(Loaders& loaders_nd,
                                         const HistAxis& ax,
                                         const SpillCut& spillcut,
                                         const Cut& cut,
                                         const Var& wei_mc,
                                         const SystShifts& shift_data,
                                         const Var& wei_data)
    : fLoadersND(loaders_nd), fAxis(ax), fSpillCut(spillcut), fCut(cut), fWeightMC(wei_mc),
      fShiftData(shift_data), fWeightData(wei_data)
  {
  }

  //----------------------------------------------------------------------
  std::unique_ptr<IPrediction> SBNExtrapGenerator::Generate(Loaders& loaders_fd,
                                                            const SystShifts& shiftMC) const
  {
    // No data shifts or weights
    return std::unique_ptr<IPrediction>(new PredictionSBNExtrap(fLoadersND, loaders_fd, fAxis, fSpillCut, fCut, shiftMC, fWeightMC, fShiftData, fWeightData));
  }
}
