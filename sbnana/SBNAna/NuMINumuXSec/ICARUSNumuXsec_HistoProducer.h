#pragma once

#include <iostream>

// ROOT
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"

// Make a plot with cuts
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

// NuMINumuXSec
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
//#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Weight.h"
//#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Systematics.h"



using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  class HistoProducer{

  public:

    HistoProducer();

    void initialize();
    bool setCut(TString cutName);

    // Spectrum booking
    // - helpers
    void FillFlashMatching(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    void FillLongestTrack(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - 221121_CRTPMTMatching
    void CRTPMTMatchingStudy(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);
    // - NuMu event selection
    void NuMINumuXSec(SpectrumLoader& loader, SpillCut spillCut=kNoSpillCut, Cut cut=kNoCut);


    void saveHistograms();

    // Systematic weights
    void setSystematicWeights();
    std::vector<const ISyst*> IAllSysts; // TODO update later
    std::vector<const ISyst*> IGENIESysts;
    std::vector< vector<Var> > vec_UniverseWeightsForEachGENIESource; // For EnsembleSpectrum
    std::vector<const ISyst*> IFluxSysts;
    std::vector<const ISyst*> IDetectorSysts;

    ~HistoProducer();

    double TargetPOT;
    TString str_TargetPOT;

    TString outputDir;
    TFile *outputfile;
    TString outputName;
    TString sampleName;

    TString currentCutName;
    vector<TString> vec_cutNames;
    std::map< TString, vector<Spectrum *> > map_cutName_to_vec_Spectrums;

    // EnsembleSpectrum-based systematics
    std::map< TString, vector< pair<TString, EnsembleSpectrum *> > > map_cutName_to_vec_SystEnsembleSpectrumPairs;
    void AddEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName);
    void AddEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName);

    // Systematics by ISyst
    std::map< TString, vector< pair<TString, Spectrum *> > > map_cutName_to_vec_SystSpectrumPairs;
    void AddUpDownSystematic(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s);
    void AddUpDownSystematic(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s);

    // booleans
    bool IsData;
    bool FillMetaData;

  private:

  };

}

