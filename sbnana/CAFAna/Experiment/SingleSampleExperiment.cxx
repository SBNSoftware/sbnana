#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"

#include "sbnana/CAFAna/Core/HistCache.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1.h"

namespace ana
{
  //----------------------------------------------------------------------
  SingleSampleExperiment::SingleSampleExperiment(const IPrediction* pred,
                                                 const Spectrum& data,
                                                 const Spectrum& cosmic)
    : fMC(pred), fData(data), fCosmic(cosmic), fMask(0)
  {
  }

  //----------------------------------------------------------------------
  SingleSampleExperiment::~SingleSampleExperiment()
  {
    delete fMask;
  }

  //----------------------------------------------------------------------
  double SingleSampleExperiment::ChiSq(osc::IOscCalcAdjustable* calc,
                                       const SystShifts& syst) const
  {
    const Spectrum p = fMC->PredictSyst(calc, syst);
    TH1D* hpred = p.ToTH1(fData.POT());

    // "Livetime" here means number of readout windows.
    if(fCosmic.Livetime() > 0){ // if cosmics supplied
      // This many are already included in the beam MC
      const double beamLivetime = p.Livetime() * fData.POT()/p.POT();
      // So we need to take this many from the cosmics-only to match the data
      const double intimeLivetime = fData.Livetime() - beamLivetime;

      TH1D* hcosmic = fCosmic.ToTH1(intimeLivetime, kLivetime);
      hpred->Add(hcosmic);
      HistCache::Delete(hcosmic);
    }

    TH1D* hdata = fData.ToTH1(fData.POT());

    // If a valid mask has been set, zero out the offending bins
    if (fMask){
      assert(hpred->GetNbinsX() == fMask->GetNbinsX());
      assert(hdata->GetNbinsX() == fMask->GetNbinsX());

      for(int i = 0; i < fMask->GetNbinsX()+2; ++i){
	if (fMask->GetBinContent(i+1) == 1) continue;
	hpred->SetBinContent(i+1, 0);
	hdata->SetBinContent(i+1, 0);
      }
    }

    const double ll = LogLikelihood(hpred, hdata);

    HistCache::Delete(hpred);
    HistCache::Delete(hdata);

    return ll;
  }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = dir;

    dir->cd();
    TObjString("SingleSampleExperiment").Write("type");

    fMC->SaveTo(dir->mkdir("mc"));
    fData.SaveTo(dir->mkdir("data"));
    fCosmic.SaveTo(dir->mkdir("cosmic"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<SingleSampleExperiment> SingleSampleExperiment::LoadFrom(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    assert(ptag->GetString() == "SingleSampleExperiment");

    assert(dir->GetDirectory("mc"));
    assert(dir->GetDirectory("data"));


    const IPrediction* mc = ana::LoadFrom<IPrediction>(dir->GetDirectory("mc")).release();
    const std::unique_ptr<Spectrum> data = Spectrum::LoadFrom(dir->GetDirectory("data"));
    const std::unique_ptr<Spectrum> cosmic = Spectrum::LoadFrom(dir->GetDirectory("cosmic"));

    auto ret = std::make_unique<SingleSampleExperiment>(mc, *data, *cosmic);
    return ret;
  }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::SetMaskHist(double xmin, double xmax, double ymin, double ymax)
  {
    fMask = GetMaskHist(fData, xmin, xmax, ymin, ymax);
  }
}
