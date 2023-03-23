#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"

HistoProducer::HistoProducer(){

  cout << "[HistoProducer::HistoProducer] called" << endl;
  TargetPOT = 6.0e20;
  str_TargetPOT = "6.0e20 POT";
  outputName = "output.root";
  currentCutName = "DefaultCutName";
  vec_cutNames.clear();
  map_cutName_to_vec_Spectrums.clear();
  map_cutName_to_vec_SystEnsembleSpectrumPairs.clear();

  FillMetaData = true;

  cout << "[HistoProducer::HistoProducer] Finished" << endl;

}

void HistoProducer::initialize(){

  cout << "[HistoProducer::initialize] outputDir = " << outputDir << endl;
  gSystem->mkdir(outputDir, kTRUE);
  outputfile = new TFile(outputDir+"/"+outputName,"RECREATE");

}

bool HistoProducer::setCut(TString cutName){

  currentCutName = cutName;
  if( std::find(vec_cutNames.begin(), vec_cutNames.end(), cutName) != vec_cutNames.end() ){
    cout << "[HistoProducer::setCut] cutName = " << cutName << " already exist.." << endl;
    return false;
  }
  else{
    cout << "[HistoProducer::setCut] adding cutName = " << cutName << endl;
    vec_cutNames.push_back(cutName);
    return true;
  }

}

void HistoProducer::FillFlashMatching(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("FMScore", Binning::Simple(102, -2., 100.), loader, varFMScore, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("FMTime", Binning::Simple(122, -62., 60.), loader, varFMTime, spillCut, cut) );

}

void HistoProducer::FillLongestTrack(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackLength", Binning::Simple(100, 0.,500.), loader, varLongestTrackLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionX", Binning::Simple(40, -1., 1.), loader, varLongestTrackDirectionX, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionY", Binning::Simple(40, -1., 1.), loader, varLongestTrackDirectionY, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackDirectionZ, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionXZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackDirectionXZ, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackForceDownDirectionX", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownDirectionX, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackForceDownDirectionY", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownDirectionY, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackForceDownDirectionZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownDirectionZ, spillCut, cut) );

}

void HistoProducer::FillSpillCount(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("CountSpill", Binning::Simple(1, 0.,1.), loader, spillvarCountSpill, spillCut)
  );

}
void HistoProducer::FillSliceCount(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("CountSlice", Binning::Simple(1, 0.,1.), loader, varCountSlice, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("IsClearCosmic", Binning::Simple(2, 0., 2.), loader, varIsClearCosmic, spillCut, cut) );

}

// - 221121_CRTPMTMatching
void HistoProducer::CRTPMTMatchingStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  FillSpillCount(loader, spillCut, cut);

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("OpFlashPeakToFirstTime", Binning::Simple(200, -100,100.), loader, spillvarOpFlashPeakToFirstTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("OpFlashTime", Binning::Simple(2500, -5.,20.), loader, ICARUSCRTPMTMatching::spillvarOpFlashTime, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("InTimeOpFlashTime", Binning::Simple(2500, -5.,20.), loader, ICARUSCRTPMTMatching::spillvarInTimeOpFlashTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TopCRTHitTime", Binning::Simple(200, -50.,50.), loader, spillvarTopCRTHitTime, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("SideCRTHitTime", Binning::Simple(200, -50.,50.), loader, spillvarSideCRTHitTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TopCRTHitPe", Binning::Simple(100, 0.,500.), loader, spillvarTopCRTHitPe, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("SideCRTHitPe", Binning::Simple(100, 0.,500.), loader, spillvarSideCRTHitPe, spillCut)
  );


  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TopCRTPMTTime", Binning::Simple(3000, -0.15, 0.15), loader, spillvarTopCRTPMTTime, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("SideCRTPMTTime", Binning::Simple(3000, -0.15, 0.15), loader, spillvarSideCRTPMTTime, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("CRTPMTTime", Binning::Simple(3000, -0.15, 0.15), loader, ICARUSCRTPMTMatching::spillvarCRTPMTTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("OpFlashTimeAllRange", Binning::Simple(6000, -3000.,3000.), loader, ICARUSCRTPMTMatching::spillvarOpFlashTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("CRTPMTMatchingID", Binning::Simple(15, -1, 14.), loader, ICARUSCRTPMTMatching::spillvarCRTPMTMatchingID, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("EastWestCRTHitPosX", Binning::Simple(300, 300.,600.), loader, spillvarEastWestCRTHitPosX, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("FlashPosX", Binning::Simple(300, 300.,600.), loader, spillvarFlashPosX, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("CRLongestTrackDirY", Binning::Simple(20, -1., 1.), loader, varCRLongestTrackDirY, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongTrackDirectionY", Binning::Simple(20, -1., 1.), loader, varLongTrackDirectionY, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionY", Binning::Simple(20, -1., 1.), loader, varLongestTrackDirectionY, spillCut, cut) );

}

// - NuMu event selection
void HistoProducer::NuMINumuXSec(SpectrumLoader& loader, SpillCut spillCut, Cut cut){
/*
  FillSpillCount(loader, spillCut, cut);
*/
  //FillFlashMatching(loader, spillCut, cut);
  FillLongestTrack(loader, spillCut, cut);

}

void HistoProducer::TruthStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  using namespace ICARUSNumuXsec::TruthMatch;

  FillSpectrum(loader, "TruthE", varNeutrinoTruthE, Binning::Simple(20, 0., 5.), spillCut, cut);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthE", Binning::Simple(20, 0., 5.), loader, varNeutrinoTruthE, spillCut, cut) );
/*
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthQ2", Binning::Simple(20, 0., 2.), loader, varTruthQ2, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Truthq0_lab", Binning::Simple(20, 0., 2.), loader, varTruthq0_lab, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Truthmodq_lab", Binning::Simple(20, 0., 2.), loader, varTruthmodq_lab, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthW", Binning::Simple(30, 0., 3.), loader, varTruthW, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthMuonKE", Binning::Simple(300, 0., 3.), loader, TruthMuonKE, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthMuonLength", Binning::Simple(500, 0., 500), loader, TruthMuonLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthProtonKE", Binning::Simple(100, 0., 1.0), loader, TruthProtonKE, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthProtonLength", Binning::Simple(100, 0., 100.), loader, TruthProtonLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionKE", Binning::Simple(100, 0., 1.0), loader, TruthChargedPionKE, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionLength", Binning::Simple(500, 0., 500.), loader, TruthChargedPionLength, spillCut, cut) );
*/
}

void HistoProducer::TrackPIDStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  using namespace ICARUSNumuXsec::TwoTrack;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTrackScore",
      Binning::Simple(20, 0., 1.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTrackScore,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackLength",
      Binning::Simple(500, 0., 500), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackDirX",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackDirX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackDirY",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackDirY,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackDirZ",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackDirZ,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2Muon",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2Muon,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2Proton",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2Proton,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackCustomChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCustomChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackCustomChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCustomChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthStartProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthStartProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthEndProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthEndProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackCollectionRR_vs_RelaxedMuonTrackCollectiondEdX", loader,
      Binning::Simple(50, 0., 25.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCollectionRR,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCollectiondEdX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthLength_vs_RelaxedMuonTrackLength", loader,
      Binning::Simple(50, 0., 500.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthLength,
      Binning::Simple(50, 0., 500.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLength,
      spillCut, cut
    )
  );

  // proton

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackLength",
      Binning::Simple(100, 0., 100), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackDirX",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackDirX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackDirY",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackDirY,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackDirZ",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackDirZ,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTrackScore",
      Binning::Simple(20, 0., 1.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTrackScore,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2Muon",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2Muon,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2Proton",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2Proton,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackCustomChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCustomChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackCustomChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCustomChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackNHitCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackNHitCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackCollectionRR_vs_RelaxedProtonTrackCollectiondEdX", loader,
      Binning::Simple(50, 0., 25.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCollectionRR,
      Binning::Simple(60, 0., 30.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCollectiondEdX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthStartProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthStartProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthEndProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthEndProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackLength_vs_RelaxedProtonTrackChi2MuonCollection", loader,
      Binning::Simple(25, 0., 50.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackLength,
      Binning::Simple(70, 0., 70.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthLength_vs_RelaxedProtonTrackLength", loader,
      Binning::Simple(100, 0., 100.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthLength,
      Binning::Simple(100, 0., 100.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackLength,
      spillCut, cut
    )
  );

}

void HistoProducer::HighChi2MuonProton(SpectrumLoader& loader, SpillCut spillCut, Cut cut){


  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackLength",
      Binning::Simple(500, 0., 500), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackLength,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2Muon",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2Muon,
      spillCut, cut
    )
  );

}

void HistoProducer::TwoTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  FillSliceCount(loader, spillCut, cut);
  FillSpillCount(loader, spillCut, cut);

  TruthStudy(loader, spillCut, cut);
  FillLongestTrack(loader, spillCut, cut);

  using namespace ICARUSNumuXsec::TwoTrack;

  TrackPIDStudy(loader, spillCut, cut);

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("CRLongestTrackDirY", Binning::Simple(100, -1., 1.), loader, varCRLongestTrackDirY, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("NPrimaryTracks", Binning::Simple(10, 0., 10.), loader, NPrimaryTracks, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackLength", Binning::Simple(500, 0., 500), loader, MuonTrackLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackP", Binning::Simple(50, 0., 5.), loader, MuonTrackP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackNuMICosineTheta", Binning::Simple(40, -1., 1.), loader, MuonTrackNuMICosineTheta, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackNuMIToVtxCosineTheta", Binning::Simple(40, -1., 1.), loader, MuonTrackNuMIToVtxCosineTheta, spillCut, cut) );
  // truth match
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonTrackTruthP_vs_MuonTrackP", loader,
      Binning::Simple(50, 0., 5.), MuonTrackTruthP,
      Binning::Simple(50, 0., 5.), MuonTrackP,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonTrackTruthNuMICosineTheta_vs_MuonTrackNuMICosineTheta", loader,
      Binning::Simple(40, -1., 1.), MuonTrackTruthNuMICosineTheta,
      Binning::Simple(40, -1., 1.), MuonTrackNuMICosineTheta,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackP", Binning::Simple(20, 0., 2.), loader, ProtonTrackP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackLength", Binning::Simple(100, 0., 100.), loader, ProtonTrackLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackNuMICosineTheta", Binning::Simple(40, -1., 1.), loader, ProtonTrackNuMICosineTheta, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackNuMIToVtxCosineTheta", Binning::Simple(40, -1., 1.), loader, ProtonTrackNuMIToVtxCosineTheta, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "ProtonTrackTruthP_vs_ProtonTrackP", loader,
      Binning::Simple(20, 0., 2.), ProtonTrackTruthP,
      Binning::Simple(20, 0., 2.), ProtonTrackP,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "ProtonTrackTruthNuMICosineTheta_vs_ProtonTrackNuMICosineTheta", loader,
      Binning::Simple(40, -1., 1.), ProtonTrackTruthNuMICosineTheta,
      Binning::Simple(40, -1., 1.), ProtonTrackNuMICosineTheta,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonProtonCosineTheta", Binning::Simple(40, -1., 1.), loader, MuonProtonCosineTheta, spillCut, cut) );

}

void HistoProducer::TwoTrackTruthMatching(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

}

void HistoProducer::ThreeTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  TwoTrackAnalysis(loader, spillCut, cut);

}

void HistoProducer::WWTPCFieldTest(SpectrumLoader& loader, SpillCut spillCut, Cut cut){


  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Z_vs_Y", loader, 
    Binning::Simple(50, 220., 420.), spillvarWWTPCTrackEndZ,
    Binning::Simple(50, 50., 150.), spillvarWWTPCTrackEndY,
    spillCut)
  );

}

void HistoProducer::Test(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackCollectionRR_vs_RelaxedMuonTrackCollectiondEdX", 
      Binning::Simple(40, 0., 40.),
      loader,
       ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCollectionRR,
      spillCut, cut
    )
  );

}

void HistoProducer::saveHistograms(){

  outputfile->cd();

  cout << "[HistoProducer::saveHistograms] Number of cuts = " << vec_cutNames.size() << endl;
  const unsigned int nCutName = vec_cutNames.size();
  double inputSamplePOT(-1.);
  double inputSampleLiveTime(-1.);
  if(nCutName>0){

    for(unsigned int ic=0; ic<nCutName; ic++){

      const TString cutName = vec_cutNames.at(ic);
      const TString dirName = cutName;
      TDirectory *dir = outputfile->GetDirectory(dirName);
      if(!dir){
        outputfile->mkdir(dirName);
        dir = outputfile->GetDirectory(dirName);
      }
      outputfile->cd(dirName);

      vector<Spectrum *> vec_Spectrums = map_cutName_to_vec_Spectrums[cutName];
      vector< pair<TString, EnsembleSpectrum *> > vec_SystEnsembleSpectrumPairs = map_cutName_to_vec_SystEnsembleSpectrumPairs[cutName];
      vector< pair<TString, Spectrum *> > vec_SystSpectrumPairs = map_cutName_to_vec_SystSpectrumPairs[cutName];

      cout << "[HistoProducer::saveHistograms] cutName = " << cutName << endl;
      cout << "[HistoProducer::saveHistograms]   Directory name = " << dirName << endl;

      // Spectrum
      cout << "[HistoProducer::saveHistograms]   Number of Spectrum = " << vec_Spectrums.size() << endl;
      for(unsigned int i=0; i<vec_Spectrums.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Spectrum.." << endl;

        if(ic==0 && i==0 && FillMetaData){
          cout << "[HistoProducer::saveHistograms]     POT = " << vec_Spectrums.at(i)->POT() << endl;
          cout << "[HistoProducer::saveHistograms]     Livetime = " << vec_Spectrums.at(i)->Livetime() << endl;
          inputSamplePOT = vec_Spectrums.at(i)->POT();
          inputSampleLiveTime = vec_Spectrums.at(i)->Livetime();
        }

        TString hName = vec_Spectrums.at(i)->GetLabels()[0];

        if(vec_Spectrums.at(i)->GetBinnings().size()==1){
          TH1 *h = vec_Spectrums.at(i)->ToTH1( vec_Spectrums.at(i)->POT() );
          cout << "[HistoProducer::saveHistograms]     Writing TH1, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }
        else if(vec_Spectrums.at(i)->GetBinnings().size()==2){
          TH2 *h = vec_Spectrums.at(i)->ToTH2( vec_Spectrums.at(i)->POT() );
          cout << "[HistoProducer::saveHistograms]     Writing TH2, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }
        else if(vec_Spectrums.at(i)->GetBinnings().size()==3){
          TH3 *h = vec_Spectrums.at(i)->ToTH3( vec_Spectrums.at(i)->POT() );
          cout << "[HistoProducer::saveHistograms]     Writing TH3, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }

      }

      // EnsembleSpectrum
      cout << "[HistoProducer::saveHistograms]   Number of SystematicEnsembleSpectrum = " << vec_SystEnsembleSpectrumPairs.size() << endl;
      for(unsigned int i=0; i<vec_SystEnsembleSpectrumPairs.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th EnsembleSpectrum.." << endl;

        TString systematicName = vec_SystEnsembleSpectrumPairs.at(i).first;
        Spectrum sNominal = vec_SystEnsembleSpectrumPairs.at(i).second->Nominal();
        TString baseLabel = sNominal.GetLabels()[0];
        dir->mkdir(baseLabel+"_"+systematicName+"_"+dirName);
        dir->cd(baseLabel+"_"+systematicName+"_"+dirName);

        if(sNominal.GetBinnings().size()==1){
          TString newLabel = baseLabel+"_NominalFromES";
          TH1 *h = sNominal.ToTH1( sNominal.POT() );
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     Nominal histogram = " << baseLabel << endl;
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH1 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH1( vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
            //TString systName = IGENIESysts.at(iu)->ShortName();
            hU->SetName(baseLabel+"_"+systematicName+"_Univ"+TString::Itoa(iu,10)+"_"+dirName);
            hU->Write();
          }
          TGraphAsymmErrors *grErrorBand = vec_SystEnsembleSpectrumPairs.at(i).second->ErrorBand(1., vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
          grErrorBand->SetName(baseLabel+"_"+systematicName+"_ErrorBandFromES_"+dirName);
          grErrorBand->Write();
        }
        else if(sNominal.GetBinnings().size()==2){
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH2 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH2( vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
            //TString systName = IGENIESysts.at(iu)->ShortName();
            hU->SetName(baseLabel+"_"+systematicName+"_Univ"+TString::Itoa(iu,10)+"_"+dirName);
            hU->Write();
          }
        }

        dir->cd();

      }

      // ISyst-based Systematic
      cout << "[HistoProducer::saveHistograms]   Number of SystematicSpectrum = " << vec_SystSpectrumPairs.size() << endl;

      for(unsigned int i=0; i<vec_SystSpectrumPairs.size(); i++){

        TString systematicName = vec_SystSpectrumPairs.at(i).first;

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Systematic Spectrum; " << systematicName << endl;

        if(vec_SystSpectrumPairs.at(i).second->GetBinnings().size()==1){
          TH1 *h = vec_SystSpectrumPairs.at(i).second->ToTH1( vec_SystSpectrumPairs.at(i).second->POT()  );
          TString baseLabel = vec_SystSpectrumPairs.at(i).second->GetLabels()[0];
          TString newLabel = baseLabel+"_"+systematicName;
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     --> Done" << endl;
        }
        else if(vec_SystSpectrumPairs.at(i).second->GetBinnings().size()==2){
          TH2 *h = vec_SystSpectrumPairs.at(i).second->ToTH2( vec_SystSpectrumPairs.at(i).second->POT() );
          TString baseLabel = vec_SystSpectrumPairs.at(i).second->GetLabels()[0];
          TString newLabel = baseLabel+"_"+systematicName;
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     --> Done" << endl;
        }

      }

      outputfile->cd();

    } // END loop cutname

    outputfile->cd();

  } // END if nCutName>0

  // When we have too many histograms, we have to run
  // Nominal and Universes separeately, and hadd later
  // For this case, some histograms should not be hadd-ed
  if(FillMetaData){

    outputfile->cd();
    outputfile->mkdir("JobInfo");
    outputfile->cd("JobInfo");
    TH1D *hist_TargetPOT = new TH1D("hist_TargetPOT", "", 1, 0., 1.);
    hist_TargetPOT->SetBinContent(1, TargetPOT);
    hist_TargetPOT->Write();

    outputfile->cd();
    outputfile->mkdir("BeamInfo");
    outputfile->cd("BeamInfo");
    TH1D *hist_POT = new TH1D("POT_BeamInfo", "POT", 1, 0., 1.);
    hist_POT->SetBinContent(1, inputSamplePOT);
    TH1D *hist_Livetime = new TH1D("Livetime_BeamInfo", "Livetime", 1, 0., 1.);
    hist_Livetime->SetBinContent(1, inputSampleLiveTime);
    hist_POT->Write();
    hist_Livetime->Write();
    outputfile->cd();
  }

}

void HistoProducer::setSystematicWeights(){

  cout << "[HistoProducer::setSystematicWeights] Setting GENIE systematics" << endl;

  const std::vector<std::string> genieKnobNames = GetSBNGenieWeightNames();
  for(const std::string& name: genieKnobNames){
    if(name=="FormZone"){
      std::cout << "[HistoProducer::setSystematicWeights] Skipping FormZone" << std::endl;
      continue;
    }
    std::string psetname = "GENIEReWeight_ICARUS_v1_multisigma_"+name;
    IGENIESysts.push_back( new SBNWeightSyst(psetname) );
  }

  cout << "[HistoProducer::setSystematicWeights] Setting flux systematics" << endl;
  IFluxSysts = GetAllNuMIFluxSysts(5);
  for(unsigned int i=0; i<IFluxSysts.size(); i++){
    cout << "[HistoProducer::setSystematicWeights] Syst = " << IFluxSysts.at(i)->ShortName() << endl;
  }

  //cout << "[HistoProducer::setSystematicWeights] Setting detector systematics" << endl;
  //IDetectorSysts.push_back( new MuonMomentumScaleSyst(0.02) );

}

void HistoProducer::AddEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName){

  map_cutName_to_vec_SystEnsembleSpectrumPairs[currentCutName].push_back( 
    std::make_pair( systName, new EnsembleSpectrum(loader, ax, spillCut, cut, varWeights) )
  );

}

void HistoProducer::AddEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName){

  map_cutName_to_vec_SystEnsembleSpectrumPairs[currentCutName].push_back(
    std::make_pair( systName, new EnsembleSpectrum(loader, axX, axY, spillCut, cut, varWeights) )
  );

}

void HistoProducer::AddUpDownSystematic(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s){

  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Up", new Spectrum(loader, ax, spillCut, cut, SystShifts(s, +1)) )
  );
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Down", new Spectrum(loader, ax, spillCut, cut, SystShifts(s, -1)) )
  );

}

void HistoProducer::AddUpDownSystematic(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s){

  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Up", new Spectrum(loader, axX, axY, spillCut, cut, SystShifts(s, +1)) )
  );
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Down", new Spectrum(loader, axX, axY, spillCut, cut, SystShifts(s, -1)) )
  );

}

void HistoProducer::FillSpectrum(SpectrumLoader& loader, const std::string& label, const Var& var, const Binning& binning, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum(label, loader, binning, var, spillCut, cut) );

  // GENIE Syst
  for(const auto& s: IGENIESysts){
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Up", new Spectrum(label, loader, binning, var, spillCut, cut, SystShifts(s, +1)) )
    );
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Down", new Spectrum(label, loader, binning, var, spillCut, cut, SystShifts(s, -1)) )
    );
  }


}

HistoProducer::~HistoProducer(){

  cout << "[HistoProducer::~HistoProducer] called" << endl;

  outputfile->Close();

  cout << "[HistoProducer::~HistoProducer] output file : " << outputDir+"/"+outputName << endl;

  cout << "[HistoProducer::~HistoProducer] Finished" << endl;

}

