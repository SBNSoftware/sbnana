#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

using namespace ana;
using namespace ICARUSNumuXsec;

void test(){

  vector<string> vec_inputs;
  for(unsigned int i=0; i<20; i++){
    TString filepath = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_68_00_01/230329_MergeMiniProd/flatcaf_"+TString::Itoa(i,10)+".root";
    vec_inputs.push_back(filepath.Data());
  }
  SpectrumLoader loader(vec_inputs);


  //SpectrumLoader loader("/pnfs//icarus/scratch/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_68_00_01/230329_MergeMiniProd/flatcaf_4.root");

  // Event (=Slilce) selection for Contained muon candidate
  Cut kCutContainedMuonPresel = kNoCut && cutRFiducial && cutNotClearCosmic && cutCRLongestTrackDirYHard && ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackContained && !ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackIsochronous && ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLengthCut;
  // Sample (based on truth-matching) selection for contained muon from nu-mu CC
  Cut kTruthContainedNuMuon = ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthContainedNuMuon;

  Spectrum* kS =  new Spectrum(
    "RelaxedMuonTrackLength", 
    Binning::Simple(500, 0., 500), 
    loader, 
    ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLength, 
    kNoSpillCut, 
    kCutContainedMuonPresel && kTruthContainedNuMuon);

  loader.Go();


  TFile *f_out = new TFile("output.root","RECREATE");
  f_out->cd();
  TH1* h = kS->ToTH1( 51770019.83189027E12 );
  h->Write();
  f_out->Close();

}
