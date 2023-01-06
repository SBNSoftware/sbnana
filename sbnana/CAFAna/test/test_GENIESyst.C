#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

using namespace ana;

//====  Var
const Var varGENIEIntCode([](const caf::SRSliceProxy* slc) -> int {
	//==== GENIE interaction code
	//==== https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360
  if(slc->truth.genie_mode<0) return -1;
  else if(slc->truth.genie_mode>13) return 14;
  else return slc->truth.genie_mode;
});
const Var varNPiBeforeFSI([](const caf::SRSliceProxy* slc) -> int {
  return slc->truth.npiminus + slc->truth.npiplus + slc->truth.npizero;
});
const Var varHitNuc([](const caf::SRSliceProxy* slc) -> int {
  return slc->truth.hitnuc;
});

//==== Cut
const Cut cutIsQE = (varGENIEIntCode==0);
const Cut cutIsRes = (varGENIEIntCode==1);
const Cut cutIsDIS = (varGENIEIntCode==2);
const Cut cutIsCC([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.iscc );
});
const Cut cutIsNC([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.isnc );
});
const Cut cutIs1Pi = (varNPiBeforeFSI==1);
const Cut cutIs2Pi = (varNPiBeforeFSI==2);
const Cut cutIsNuMu([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.pdg==14 );
});
const Cut cutIsNuMuBar([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.pdg==-14 );
});
const Cut cutHitProton = kIsNuSlice && (varHitNuc==2212);
const Cut cutHitNeutron = kIsNuSlice && (varHitNuc==2112);

//==== fwd decl
void setCanvas(TCanvas* c1, TPad* c1_up, TPad* c1_down);
void setAxis(TH1D* hist_up, TH1D* hist_down);

void test_GENIESyst(){

  //==== Input
  SpectrumLoader loader("/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_37_01_03p02/MultiSigmaAdded_RenamedPSet/*.root");

  //==== output; where plots are saved
  TString outputPath = "./test_GENIESyst/";
  gSystem->mkdir(outputPath, kTRUE);

  //==== Variable
  const Var kTrueE = SIMPLEVAR(truth.E);

  //==== (Truth) Cut

  std::vector<Cut> vec_Cuts = {
    kIsNuSlice, 
    cutIsQE,
    cutIsRes,
    cutIsDIS,
    //==== below for NonResBG
    cutIsDIS && cutIsNuMuBar && cutHitNeutron && cutIsCC,
    cutIsDIS && cutIsNuMuBar && cutHitNeutron && cutIsNC,
    cutIsDIS && cutIsNuMuBar && cutHitProton && cutIsCC,
    cutIsDIS && cutIsNuMuBar && cutHitProton && cutIsNC,
    cutIsDIS && cutIsNuMu && cutHitNeutron && cutIsCC,
    cutIsDIS && cutIsNuMu && cutHitNeutron && cutIsNC,
    cutIsDIS && cutIsNuMu && cutHitProton && cutIsCC,
    cutIsDIS && cutIsNuMu && cutHitProton && cutIsNC,
  };
  std::vector<std::string> vec_CutNames = {
    "All",
    "QE",
    "Res",
    "DIS",
    //==== below for NonResBG
    "DIS_vbarnCC",
    "DIS_vbarnNC",
    "DIS_vbarpCC",
    "DIS_vbarpNC",
    "DIS_vnCC",
    "DIS_vnNC",
    "DIS_vpCC",
    "DIS_vpNC",
  };
  std::vector<std::string> vec_CutLatexNames = {
    "All interactions",
    "QE",
    "Resonance",
    "DIS",
    //==== below for NonResBG
    "DIS, #bar{#nu_{#mu}}, n, CC",
    "DIS, #bar{#nu_{#mu}}, n, NC",
    "DIS, #bar{#nu_{#mu}}, p, CC",
    "DIS, #bar{#nu_{#mu}}, p, NC",
    "DIS, #nu_{#mu}, n, CC",
    "DIS, #nu_{#mu}, n, NC",
    "DIS, #nu_{#mu}, p, CC",
    "DIS, #nu_{#mu}, p, NC",
  };

  //==== Binning
  const Binning binsEnergy = Binning::Simple(30, 0, 6);
  const HistAxis axEnergy("True neutrino energy (GeV)", binsEnergy, kTrueE);

  //==== Spectrum maps
  std::map<TString, Spectrum*> map_Spectrums;
  std::map<TString, EnsembleSpectrum*> map_EnsembleSpectrums;

  //==== Systematics; multisigma
  const std::vector<const ISyst*> IGENIEMultisigmaSysts = GetSBNGenieWeightSysts();
  std::vector<Var> weis;
  weis.reserve(1000);
  for(int i = 0; i < 1000; ++i) weis.push_back(GetUniverseWeight("multisim_Genie", i));

  //==== Loop over cuts and define (Ensemble)Spectrums
  for(unsigned int i_cut=0; i_cut<vec_CutNames.size(); i_cut++){

    Cut this_Cut = vec_Cuts.at(i_cut);
    std::string this_CutName = vec_CutNames.at(i_cut);

    //==== CV
    map_Spectrums[this_CutName+"_CV"] = new Spectrum(loader, axEnergy, kNoSpillCut, this_Cut);

    //==== Multisigma
    for(unsigned int i_syst=0; i_syst<IGENIEMultisigmaSysts.size(); i_syst++){
      const ISyst* iS = IGENIEMultisigmaSysts.at(i_syst);
      map_Spectrums[this_CutName+"_"+iS->LatexName()+"_1Up"] = new Spectrum(loader, axEnergy, kNoSpillCut, this_Cut, SystShifts(iS, +1));
      map_Spectrums[this_CutName+"_"+iS->LatexName()+"_1Dn"] = new Spectrum(loader, axEnergy, kNoSpillCut, this_Cut, SystShifts(iS, -1));
      map_Spectrums[this_CutName+"_"+iS->LatexName()+"_2Up"] = new Spectrum(loader, axEnergy, kNoSpillCut, this_Cut, SystShifts(iS, +2));
      map_Spectrums[this_CutName+"_"+iS->LatexName()+"_2Dn"] = new Spectrum(loader, axEnergy, kNoSpillCut, this_Cut, SystShifts(iS, -2));
      map_Spectrums[this_CutName+"_"+iS->LatexName()+"_3Up"] = new Spectrum(loader, axEnergy, kNoSpillCut, this_Cut, SystShifts(iS, +3));
      map_Spectrums[this_CutName+"_"+iS->LatexName()+"_3Dn"] = new Spectrum(loader, axEnergy, kNoSpillCut, this_Cut, SystShifts(iS, -3));
    }

    //==== Multisim; use EnsembleSpectrum
    map_EnsembleSpectrums[this_CutName+"_multisim"] = new EnsembleSpectrum(loader, axEnergy, kNoSpillCut, this_Cut, weis);

  }

  //==== Go
  loader.Go();

  //===============
  //==== Plotting
  //===============

  //==== y=1 line for the ratio
  double x_1[2], y_1[2];
  x_1[0] = 9000;  y_1[0] = 1;
  x_1[1] = -9000;  y_1[1] = 1;
  TGraph *g1 = new TGraph(2, x_1, y_1);

  for(unsigned int i_cut=0; i_cut<vec_Cuts.size(); i_cut++){

    TString this_CutName = vec_CutNames.at(i_cut);
    TString this_CutLatexName = vec_CutLatexNames.at(i_cut);

    //==== CV
    TH1D *h_CV = map_Spectrums[this_CutName+"_CV"]->ToTH1(kPOTnominal);
    h_CV->SetLineColor(kBlack);

    //==== TLatex
    TLatex tl;
    tl.SetNDC();
    tl.SetTextSize(0.035);

    //==== Multisigma
    for(unsigned int i_syst=0; i_syst<IGENIEMultisigmaSysts.size(); i_syst++){

      TString knobName = IGENIEMultisigmaSysts.at(i_syst)->LatexName();

      TH1D *h_1Up = map_Spectrums[this_CutName+"_"+knobName+"_1Up"]->ToTH1(kPOTnominal);
      TH1D *h_1Dn = map_Spectrums[this_CutName+"_"+knobName+"_1Dn"]->ToTH1(kPOTnominal);
      TH1D *h_2Up = map_Spectrums[this_CutName+"_"+knobName+"_2Up"]->ToTH1(kPOTnominal);
      TH1D *h_2Dn = map_Spectrums[this_CutName+"_"+knobName+"_2Dn"]->ToTH1(kPOTnominal);
      TH1D *h_3Up = map_Spectrums[this_CutName+"_"+knobName+"_3Up"]->ToTH1(kPOTnominal);
      TH1D *h_3Dn = map_Spectrums[this_CutName+"_"+knobName+"_3Dn"]->ToTH1(kPOTnominal);

      h_1Up->SetLineColor(kBlue);
      h_1Dn->SetLineColor(kBlue);
      h_2Up->SetLineColor(kRed);
      h_2Dn->SetLineColor(kRed);
      h_3Up->SetLineColor(kViolet);
      h_3Dn->SetLineColor(kViolet);

      h_1Dn->SetLineStyle(2);
      h_2Dn->SetLineStyle(2);
      h_3Dn->SetLineStyle(2);

      //==== Legend
      TLegend *lg = new TLegend(0.62, 0.62, 0.95, 0.95);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);
      lg->AddEntry(h_CV, "CV", "l");

      lg->AddEntry(h_1Up, knobName+", +1#sigma", "l");
      lg->AddEntry(h_1Dn, knobName+", -1#sigma", "l");
      lg->AddEntry(h_2Up, knobName+", +2#sigma", "l");
      lg->AddEntry(h_2Dn, knobName+", -2#sigma", "l");
      lg->AddEntry(h_3Up, knobName+", +3#sigma", "l");
      lg->AddEntry(h_3Dn, knobName+", -3#sigma", "l");

      TCanvas *c = new TCanvas("c_"+this_CutName+"_"+knobName, "", 800, 800);
      TPad *c_rate = new TPad("c_rate_"+this_CutName+"_"+knobName, "", 0, 0.25, 1, 1);
      TPad *c_ratio = new TPad("c_ratio_"+this_CutName+"_"+knobName, "", 0, 0, 1, 0.25);

      setCanvas(c, c_rate, c_ratio);

      c->cd();
      c_rate->Draw();
      c_ratio->Draw();

      //==== rate
      c_rate->cd();

      h_1Up->Draw("histsame");
      h_1Dn->Draw("histsame");
      h_2Up->Draw("histsame");
      h_2Dn->Draw("histsame");
      h_3Up->Draw("histsame");
      h_3Dn->Draw("histsame");
      h_CV->Draw("histsame");

      h_1Up->GetYaxis()->SetRangeUser(0., 1.1*max(h_3Up->GetMaximum(), h_3Dn->GetMaximum()));

      //==== ratio
      c_ratio->cd();

      TH1D *h_1UpRatio = (TH1D *)h_1Up->Clone("h_1UpRatio");
      TH1D *h_1DnRatio = (TH1D *)h_1Dn->Clone("h_1DnRatio");
      TH1D *h_2UpRatio = (TH1D *)h_2Up->Clone("h_2UpRatio");
      TH1D *h_2DnRatio = (TH1D *)h_2Dn->Clone("h_2DnRatio");
      TH1D *h_3UpRatio = (TH1D *)h_3Up->Clone("h_3UpRatio");
      TH1D *h_3DnRatio = (TH1D *)h_3Dn->Clone("h_3DnRatio");

      h_1UpRatio->Divide(h_CV);
      h_1DnRatio->Divide(h_CV);
      h_2UpRatio->Divide(h_CV);
      h_2DnRatio->Divide(h_CV);
      h_3UpRatio->Divide(h_CV);
      h_3DnRatio->Divide(h_CV);

      h_1UpRatio->Draw("histsame");
      h_1DnRatio->Draw("histsame");
      h_2UpRatio->Draw("histsame");
      h_2DnRatio->Draw("histsame");
      h_3UpRatio->Draw("histsame");
      h_3DnRatio->Draw("histsame");

      h_1UpRatio->GetYaxis()->SetTitle("Ratio to CV");
      h_1UpRatio->GetYaxis()->SetRangeUser(0.5, 1.5);

      g1->Draw("same");

      //==== Finalize
      c->cd();

      setAxis(h_1Up, h_1UpRatio);

      lg->Draw();
      tl.DrawLatex(0.30, 0.96, TString::Format("%1.2e POT",kPOTnominal) );
      tl.DrawLatex(0.60, 0.96, "Interaction : "+this_CutLatexName);

      c->SaveAs(outputPath+this_CutName+"_"+knobName+".pdf");
      c->Close();

    } // END Multisigma loop

    //==== Multisim

    EnsembleSpectrum* es = map_EnsembleSpectrums[this_CutName+"_multisim"];

    TCanvas *c_es = new TCanvas("c_es_"+this_CutName, "", 800, 800);
    TPad *c_es_rate = new TPad("c_es_rate_"+this_CutName, "", 0, 0.25, 1, 1);
    TPad *c_es_ratio = new TPad("c_es_ratio_"+this_CutName, "", 0, 0, 1, 0.25);

    setCanvas(c_es, c_es_rate, c_es_ratio);

    c_es->cd();
    c_es_rate->Draw();
    c_es_ratio->Draw();

    //==== Legend
    TLegend *lg_es = new TLegend(0.67, 0.80, 0.92, 0.92);
    lg_es->SetBorderSize(0);
    lg_es->SetFillStyle(0);
    lg_es->AddEntry(h_CV, "CV", "l");

    //==== rate

    c_es_rate->cd();
    h_CV->Draw("axis");

    double y_max(-99);
    for(unsigned int i = 0; i < es->NUniverses(); ++i){
      TH1D *h_es = es->Universe(i).ToTH1(kPOTnominal);
      h_es->SetLineColor(kRed-10);
      if(i==0) lg_es->AddEntry(h_es, "Universes", "l");
      y_max = max(y_max, h_es->GetMaximum());

      c_es_rate->cd();
      h_es->Draw("histsame");

    }

    TGraphAsymmErrors* g_errband = es->ErrorBand(kPOTnominal, +1.);
    g_errband->Draw("esame");

    h_CV->Draw("histsame");
    h_CV->GetYaxis()->SetRangeUser(0., 1.1*y_max);

    //==== ratio

    c_es_ratio->cd();
    TH1D *h_axis_down = new TH1D("h_axis_down_"+this_CutName, "", binsEnergy.NBins(), binsEnergy.Min(), binsEnergy.Max());
    h_axis_down->GetYaxis()->SetRangeUser(0.5,1.5);
    h_axis_down->SetTitle("");
    h_axis_down->GetYaxis()->SetTitle("Ratio to CV");
    h_axis_down->GetXaxis()->SetTitle("True neutrino energy (GeV)");
    h_axis_down->Draw("axis");

    TGraphAsymmErrors* g_errband_ratio = (TGraphAsymmErrors*)g_errband->Clone("g_errband_ratio");

    for(int i=0; i<g_errband_ratio->GetN(); i++){
      double x, y, yerr_low, yerr_high;

      g_errband_ratio->GetPoint(i, x, y);
      yerr_low  = g_errband_ratio->GetErrorYlow(i);
      yerr_high = g_errband_ratio->GetErrorYhigh(i);

      g_errband_ratio->SetPoint(i, x, 1);
      g_errband_ratio->SetPointEYlow(i, yerr_low/y);
      g_errband_ratio->SetPointEYhigh(i, yerr_high/y);
    }
    g_errband_ratio->Draw("psame");
    g1->Draw("same");

    h_axis_down->GetYaxis()->SetRangeUser(0.5,1.5);

    //==== finalize
    c_es->cd();

    setAxis(h_CV, h_axis_down);

    lg_es->Draw();

    tl.DrawLatex(0.30, 0.96, TString::Format("%1.2e POT",kPOTnominal) );
    tl.DrawLatex(0.60, 0.96, "Interaction : "+this_CutLatexName);

    c_es->SaveAs(outputPath+this_CutName+"_Multisim.pdf");
    c_es->Close();

  } // END cut loop

}


void setCanvas(TCanvas *c1, TPad *c1_up, TPad *c1_down){

  c1->SetTopMargin( 0.05 );
  c1->SetBottomMargin( 0.13 );
  c1->SetRightMargin( 0.05 );
  c1->SetLeftMargin( 0.16 );

  c1_up->SetTopMargin( 0.07 );
  c1_up->SetBottomMargin( 0.025 );
  c1_up->SetLeftMargin( 0.15 );
  c1_up->SetRightMargin( 0.032 );

  c1_down->SetTopMargin( 0.035 );
  c1_down->SetBottomMargin( 0.4 );
  c1_down->SetLeftMargin( 0.15 );
  c1_down->SetRightMargin( 0.032 );

}

void setAxis(TH1D* hist_up, TH1D* hist_down){

  hist_up->GetYaxis()->SetLabelSize(0.05);
  hist_up->GetYaxis()->SetTitleSize(0.070);
  hist_up->GetYaxis()->SetTitleOffset(1.10);
  hist_up->GetXaxis()->SetLabelSize(0);

  hist_down->SetTitle("");

  hist_down->GetXaxis()->SetLabelSize(0.13);
  hist_down->GetXaxis()->SetTitleSize(0.18);

  hist_down->GetYaxis()->SetLabelSize(0.13);
  hist_down->GetYaxis()->SetTitleSize(0.16);
  hist_down->GetYaxis()->SetTitleOffset(0.4);

}



