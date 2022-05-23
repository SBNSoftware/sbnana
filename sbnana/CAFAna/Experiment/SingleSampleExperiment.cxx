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
                                                 const Spectrum& cosmicInTime,
                                                 const Spectrum& cosmicOutOfTime)
    : fMC(pred), fData(data), fCosmicInTime(cosmicInTime), fCosmicOutOfTime(cosmicOutOfTime), fMask(0)
  {
    if(cosmicInTime.POT() > 0){
      std::cout << "SingleSampleExperiment: in-time cosmics have nonzero POT. "
                << "Did you confuse the two cosmics arguments?"
                << std::endl;
      abort();
    }
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
    Spectrum p = fMC->PredictSyst(calc, syst);

    if(fCosmicOutOfTime.POT() > 0){ // if out-of-time cosmics supplied
      p += fCosmicOutOfTime;
    }

    TH1D* hpred = p.ToTH1(fData.POT());

    // "Livetime" here means number of readout windows.
    if(fCosmicInTime.Livetime() > 0){ // if in-cosmics supplied
      // This many are already included in the beam MC
      const double beamLivetime = p.Livetime() * fData.POT()/p.POT();
      // So we need to take this many from the cosmics-only to match the data
      const double intimeLivetime = fData.Livetime() - beamLivetime;

      TH1D* hcosmic = fCosmicInTime.ToTH1(intimeLivetime, kLivetime);
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
    fCosmicInTime.SaveTo(dir->mkdir("cosmicInTime"));
    fCosmicOutOfTime.SaveTo(dir->mkdir("cosmicOutOfTime"));

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

    // Legacy format. This is the in-time cosmics
    if(dir->GetDirectory("cosmic")){
      const std::unique_ptr<Spectrum> cosmicInTime = Spectrum::LoadFrom(dir->GetDirectory("cosmic"));
      return std::make_unique<SingleSampleExperiment>(mc, *data, *cosmicInTime, Spectrum(0, {}, {}, 0, 0));
    }

    const std::unique_ptr<Spectrum> cosmicInTime = Spectrum::LoadFrom(dir->GetDirectory("cosmicInTime"));
    const std::unique_ptr<Spectrum> cosmicOutOfTime = Spectrum::LoadFrom(dir->GetDirectory("cosmicOutOfTime"));

    return std::make_unique<SingleSampleExperiment>(mc, *data, *cosmicInTime, *cosmicOutOfTime);
  }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::SetMaskHist(double xmin, double xmax, double ymin, double ymax)
  {
    fMask = GetMaskHist(fData, xmin, xmax, ymin, ymax);
  }
}
