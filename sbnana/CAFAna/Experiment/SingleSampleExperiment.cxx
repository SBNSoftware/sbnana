#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"

#include "sbnana/CAFAna/Core/LoadFromFile.h"
//#include "CAFAna/Core/StanUtils.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH2.h"

namespace ana
{
  //----------------------------------------------------------------------
  SingleSampleExperiment::SingleSampleExperiment(const IPrediction* pred,
                                                 const Spectrum& data,
                                                 const Spectrum& cosmic)
    : fMC(pred), fData(data), fCosmic(cosmic)
  {
  }

  //----------------------------------------------------------------------
  SingleSampleExperiment::~SingleSampleExperiment()
  {
  }

  //----------------------------------------------------------------------
  double SingleSampleExperiment::ChiSq(osc::IOscCalcAdjustable* calc,
                                       const SystShifts& syst) const
  {
    const Spectrum pred = fMC->PredictSyst(calc, syst);
    Eigen::ArrayXd apred = pred.GetEigen(fData.POT());
    Eigen::ArrayXd adata = fData.GetEigen(fData.POT());

    // "Livetime" here means number of readout windows.
    if(fCosmic.Livetime() > 0){ // if cosmics supplied
      // This many are already included in the beam MC
      const double beamLivetime = pred.Livetime() * fData.POT()/pred.POT();
      // So we need to take this many from the cosmics-only to match the data
      const double intimeLivetime = fData.Livetime() - beamLivetime;

      Eigen::ArrayXd acosmic = fCosmic.GetEigen(intimeLivetime, kLivetime);
      apred += acosmic;
    }

    ApplyMask(apred, adata);

    // full namespace qualification to avoid degeneracy with method inherited
    // from IExperiment
    return ana::LogLikelihood(apred, adata);
   }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::ApplyMask(Eigen::ArrayXd& a,
                                         Eigen::ArrayXd& b) const
  {
    if(fMaskA.size() == 0) return;

    assert(a.size() == fMaskA.size());
    assert(b.size() == fMaskA.size());

    // Arrays mean we get bin-by-bin operations
    a *= fMaskA;
    b *= fMaskA;
  }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = dir;

    dir->cd();
    TObjString("SingleSampleExperiment").Write("type");

    fMC->SaveTo(dir->mkdir("mc"));
    fData.SaveTo(dir, "data");
    fCosmic.SaveTo(dir, "cosmic");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<SingleSampleExperiment> SingleSampleExperiment::LoadFrom(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    assert(ptag->GetString() == "SingleSampleExperiment");

    assert(dir->GetDirectory("mc"));


    const IPrediction* mc = ana::LoadFrom<IPrediction>(dir->GetDirectory("mc")).release();
    const std::unique_ptr<Spectrum> data = Spectrum::LoadFrom(dir, "data");
    const std::unique_ptr<Spectrum> cosmic = Spectrum::LoadFrom(dir, "cosmic");

    auto ret = std::make_unique<SingleSampleExperiment>(mc, *data);

    return ret;
  }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::SetMaskHist(double xmin, double xmax, double ymin, double ymax)
  {
    fMaskA = GetMaskArray(fData, xmin, xmax, ymin, ymax);
  }
}
