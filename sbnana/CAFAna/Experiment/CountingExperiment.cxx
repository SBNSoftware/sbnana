#include "sbnana/CAFAna/Experiment/CountingExperiment.h"

#include "sbnana/CAFAna/Prediction/IPrediction.h"

#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //----------------------------------------------------------------------
  CountingExperiment::~CountingExperiment()
  {
  }

  //----------------------------------------------------------------------
  double CountingExperiment::ChiSq(osc::IOscCalcAdjustable* osc,
                                   const SystShifts& syst) const
  {
    double exp = fMC->PredictSyst(osc, syst).Integral(fData.POT());
    double obs = fData.Integral(fData.POT());

    return LogLikelihood(exp, obs);
  }

  //----------------------------------------------------------------------
  void CountingExperiment::SaveTo(TDirectory* dir, const std::string& name) const
  {
    TDirectory* tmp = dir;

    dir = dir->mkdir(name.c_str()); // switch to subdir
    dir->cd();

    TObjString("CountingExperiment").Write("type");

    fMC->SaveTo(dir, "mc");
    fData.SaveTo(dir, "data");

    dir->Write();
    delete dir;

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<CountingExperiment> CountingExperiment::LoadFrom(TDirectory* dir, const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str()); // switch to subdir
    assert(dir);

    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    assert(ptag->GetString() == "CountingExperiment");

    assert(dir->GetDirectory("mc"));
    

    const IPrediction* mc = ana::LoadFrom<IPrediction>(dir, "mc").release();
    const std::unique_ptr<Spectrum> data = Spectrum::LoadFrom(dir, "data");

    auto ret = std::make_unique<CountingExperiment>(mc, *data);
    return ret;
  }
}
