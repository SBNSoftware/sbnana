// Get Numbers for various interaction types and flavours
// cafe event_numbers.C

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

// #include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"
#include "sbnana/CAFAna/Analysis/Calcs.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "OscLib/OscCalcSterile.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"

using namespace ana;

// Put an "SBND simulation" tag in the corner                                                       
void Experiment(std::string expt)
{
  TLatex* prelim = new TLatex(.9, .95, (expt+" Simulation").c_str());
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

void event_numbers_finalstate(const std::string expt = "SBND")
{
  const std::string dir = "/sbnd/data/users/jlarkin/workshop_samples/";
  const std::string fnameBeam = dir + "output_SBNOsc_NumuSelection_Modern_" + expt + ".flat.root";

  // Source of events
  SpectrumLoader loaderBeam(fnameBeam);
  // SpectrumLoader loaderSwap(fnameSwap);

  // Vars and cuts

  // a struct to hold the different flavours and signs
  const unsigned int kNumFlavours = 3;
  const unsigned int kNumSigns    = 2;
  const unsigned int kNumCurrents = 2;

  struct Flavour
  {
    ana::Flavors::Flavors_t flav;
    string label;
  };
  Flavour flav[kNumFlavours] = {
    {ana::Flavors::kAllNuMu, "numu"},
    {ana::Flavors::kAllNuE,  "nue"},
    {ana::Flavors::kAll,  "all"}
  };
  struct Signtype
  {
    ana::Sign::Sign_t sign;
    string label;
  };
    Signtype sign[kNumSigns] = {
    {ana::Sign::kNu,     ""},
    {ana::Sign::kAntiNu, "anti"}
  };
  struct Currtype
  {
    ana::Current::Current_t curr;
    string label;
  };
    Currtype curr[kNumCurrents] = {
    {ana::Current::kCC, "CC"},
    {ana::Current::kNC, "NC"}
  };

  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  const Var kWeight = SIMPLEVAR(reco.weight);

  const Var kNumPiZeros(
         [](const caf::SRSliceProxy *sr)
         {
          unsigned int nPis = 0;
          for(unsigned int part_idx = 0; part_idx < sr->truth[0].finalstate.size(); part_idx++){
            if(sr->truth[0].finalstate[part_idx].pdg == 111) nPis++;
           }
          return nPis;
          });

  const Var kNumPiCharged(
           [](const caf::SRSliceProxy *sr)
           {
            unsigned int nPis = 0;
            for(unsigned int part_idx = 0; part_idx < sr->truth[0].finalstate.size(); part_idx++){
              if( abs(sr->truth[0].finalstate[part_idx].pdg) == 211) nPis++;
             }
            return nPis;
            });

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  const Var kCC = SIMPLEVAR(truth[0].neutrino.iscc);
  const Cut kIsCC = kCC > 0;

  const Cut kZeroPi  = ((kNumPiZeros + kNumPiCharged) == 0);
  const Cut kOnePiCh = (kNumPiZeros == 0 && kNumPiCharged == 1);
  const Cut kOnePi0  = (kNumPiZeros == 1 && kNumPiCharged == 0);
  const Cut kTwoPiCh = (kNumPiZeros == 0 && kNumPiCharged == 2);
  const Cut kOneEach = (kNumPiZeros == 1 && kNumPiCharged == 1);
  const Cut kTwoPi0  = (kNumPiZeros == 2 && kNumPiCharged == 0);
  const Cut kElsePi  = (!kZeroPi && !kOnePiCh && !kOnePi0 
                        && !kTwoPiCh && !kOneEach && !kTwoPi0);

  const unsigned int kNumModes = 8;
  struct Selection
  {
    Cut cut;
    string label;
    int colour;
  };
  Selection cuts[kNumModes] = {
    {kZeroPi,  "0pi",   kRed},
    {kOnePiCh, "1piCh", kGreen+2},
    {kOnePi0,  "1pi0",  kBlue},
    {kTwoPiCh, "2piCh", kViolet},
    {kOneEach, "1each", kGray+2},
    {kTwoPi0,  "2pi0",  kCyan+2},
    {kElsePi,  "else",  kOrange+3},
    {kNoCut,   "All",   kBlack}
  };

  Spectrum kNumPiZerosSpec("Number of pions", Binning::Simple(5,0,5), loaderBeam, kNumPiZeros, kNoCut, kNoShift, kWeight);
  Spectrum kNumPiChargedSpec("Number of pions", Binning::Simple(5,0,5), loaderBeam, kNumPiCharged, kNoCut, kNoShift, kWeight);

  // For now let's not bother with oscillations
  osc::NoOscillations* noosc = new osc::NoOscillations;

  // The main prediction
  PredictionNoExtrap* pred[kNumModes];

  for(unsigned int m = 0; m < kNumModes; m++){
    pred[m] = new PredictionNoExtrap(loaderBeam, kNullLoader, kNullLoader, kNullLoader,
      axEnergy, cuts[m].cut, kNoShift, kWeight);
  }

  // GO!
  loaderBeam.Go();
  // loaderSwap.Go();

  const double pot = kPOTnominal;

  std::cout << "     NuMu,  Nue,  All,  anti-Numu,  anti-Nue,  anti-All,  NC" << std::endl;

  TH1* h[kNumModes][kNumCurrents][kNumSigns][kNumFlavours];

  for(unsigned int m = 0; m < kNumModes; m++){

    std::cout << cuts[m].label;

    for(unsigned int k = 0; k < kNumCurrents; k++){
      for(unsigned int j = k; j < kNumSigns; j++){ // loop over signs only for CC
        for(unsigned int i = k; i < kNumFlavours-k; i++){ // loop over flavours only for CC
          h[m][k][j][i] = pred[m]->PredictComponent(noosc, k == 0 ? flav[i].flav : Flavors::kAll, 
                                    curr[k].curr, (k == 0) ? sign[j].sign : Sign::kBoth).ToTH1(pot);
          float nEvts = h[m][k][j][i]->Integral();
          std::cout << ", " << nEvts ;
        } // end loop over flavours
      } // end loop over signs
    } // end loop over currents
    std::cout << endl;
  } // end loop over modes

  TFile* fOutput = new TFile(("output/output_finalstate_"+expt+".root").c_str(),"RECREATE");

  TCanvas* c1 = new TCanvas("c1","c1");

  for(unsigned int m = 0; m < kNumModes; m++){
    for(unsigned int k = 0; k < kNumCurrents; k++){ 
      for(unsigned int j = k; j < kNumSigns; j++){ // loop over signs only for CC
        for(unsigned int i = k; i < kNumFlavours-k; i++){ // loop over flavours only for CC

          string currLabel = curr[k].label;
          string flavLabel = ( (k == 0 && i != 2) ? flav[i].label : "all");
          string signLabel = ( (k == 0) ? sign[j].label : "");

          h[m][k][j][i]->SetLineColor(cuts[m].colour);
          h[m][k][j][i]->SetMarkerColor(cuts[m].colour);
          h[m][k][j][i]->Write((expt+"_"+cuts[m].label+"_"+currLabel+"_"+signLabel+flavLabel).c_str());   
          h[m][k][j][i]->SetName((cuts[m].label+"_"+currLabel+"_"+signLabel+flavLabel).c_str());      
          h[m][k][j][i]->Draw("hist");
          Experiment(expt);
          c1->SaveAs(("output/"+expt+"_"+cuts[m].label+"_"+currLabel+"_"+signLabel+flavLabel+".pdf").c_str());

        } // loop over flavours
      } // loop over signs
    } // loop over currents
  } // loop over modes


  TH1* hNumPiZeros = kNumPiZerosSpec.ToTH1(pot);
  hNumPiZeros->SetLineColor(kRed);
  hNumPiZeros->SetFillColor(kRed);
  TH1* hNumPiCharged = kNumPiChargedSpec.ToTH1(pot);
  hNumPiCharged->SetLineColor(kBlue);
  hNumPiCharged->SetFillColor(kBlue);

  hNumPiZeros->Write("hNumPiZeros");
  hNumPiCharged->Write("hNumPiCharged");

  fOutput->Close();

}
