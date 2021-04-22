#include "CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"

#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TText.h"

#include <iostream>
#include <fstream>
#include <iomanip>

ofstream output("/outputdir_path/table_output.txt");

namespace ana{

  // ----------------------------------------------------------------------
  // Tables

  // Highlight cell depending on range
  TString thisCellColor(double weird){
    TString mystring = "";
    double thisval = abs(weird);
    if( thisval>=0.05 && thisval<0.10) mystring = "\\cellcolor{green!25}";
    else if( thisval>=0.10 && thisval<0.15) mystring = "\\cellcolor{yellow!25}";
    else if( thisval>=0.15) mystring = "\\cellcolor{red!25}";
    return mystring;
  }

  TString fixLatexName(TString mystring){
    // latex doesnt like underscores
    std::vector<TString> in = {"#"," ",".","_"};
    std::vector<TString> out = {"","","","\\_"};

    for(unsigned int i=0;i<in.size();i++)
      mystring.ReplaceAll(in[i],out[i]);
    return mystring;
  }

  void printTableHeader(int quantId=0)
  {
    std::setprecision(3);
    output << "\\begin{table}[H]\n";
    output << "\\centering\n";
    output << "\\resizebox{\\textwidth}{!}{\n";
    // output << "\\begin{tabular}{|l||c|c|c|c|c||c|c|}\n";
    output << "\\begin{tabular}{|l||c|c|c|c|c||c|}\n";
    output << "\\hline \n";
    // output << "\\multicolumn{1}{|c||}{} & \\multicolumn{5}{c||}{Number of interactions (\\% total)} & \\multicolumn{2}{c|}{Integrated} \\\\ \\hline \n \
    // \\multicolumn{1}{|c||}{Cut} & $\\nu_{e}$ CC & $\\nu_{\\mu}$ CC & NC & Cosmic & Other bkg & Efficiency & Purity \\\\ \\hline \n";
    output << "\\multicolumn{1}{|c||}{} & \\multicolumn{4}{c||}{Number of interactions (\\% total)} & \\multicolumn{2}{c|}{Integrated} \\\\ \\hline \n \
    \\multicolumn{1}{|c||}{Cut} & $\\nu_{e}$ CC & $\\nu_{\\mu}$ CC & NC & Cosmic & Efficiency & Purity \\\\ \\hline \n";
  }// printTableHeader

  void printTableFooter(){
    output << "\\end{tabular}}\n";
    output << "\\end{table}";
    output << "\n\n\n";
    //output.close();
  }

  void printEventsLine(std::string cutname, float nue, float numu, float nc, float cos, float other, float eff, float pur){
    float total = nue+numu+nc+cos+other;
    float percnue    = 100*nue/total;
    float percnumu   = 100*numu/total;
    float percnc     = 100*nc/total;
    float perccos    = 100*cos/total;
    float percother  = 100*other/total;
    float perctotbkg = 100*(numu+nc+cos+other)/total;
    // output << std::fixed << std::setw(6) << std::setprecision(3) << cutname << "&" << nue <<"("<<percnue<<")"<< "&" << numu<<"("<<percnumu<<")"<< "&" << nc<<"("<<percnc<<")"<< "&" << cos<<"("<<perccos<<")"<< "&" << other<<"("<<percother<<")"<< "&" << eff << "&" << pur << "\\\\ \\hline \n";
    output << std::fixed << std::setw(6) << std::setprecision(3) << cutname << "&" << nue <<"("<<percnue<<")"<< "&" << numu<<"("<<percnumu<<")"<< "&" << nc<<"("<<percnc<<")"<< "&" << cos<<"("<<perccos<<")"<< "&" << eff << "&" << pur << "\\\\ \\hline \n";
  }


  // ----------------------------------------------------------------------
  // Canvases
  // split canvas in 2
  void SplitCanvas2(TCanvas *& c1, TPad *& pad1, TPad *& pad2){

    c1 = new TCanvas("c1","",700,800);
    c1->cd();

    pad1 = new TPad("pad1","pad1",0,0,1,1);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0.4);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.03);
    pad1->SetFillStyle(0);
    pad1->Draw();
    c1->cd();

    pad2 = new TPad("pad2","pad2",0,0,1,1);// x1 y1 x2 y2
    pad2->SetTopMargin(0.6);
    pad2->SetBottomMargin(0.1);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.03);
    pad2->SetFillStyle(0);
    pad2->Draw();
    c1->cd();

  }

  float GetHistMax(std::vector<TH1*> histos){

    float hmax = 0.;
    for(unsigned int hId=0; hId<histos.size(); hId++){
      float thismax = histos[hId]->GetMaximum();
      if(thismax>hmax) hmax=thismax;
    }
    return hmax;
  }


  void PimpHist(TH1* histo, Color_t color, Style_t linestyle, int linewidth, Style_t markerstyle=8, double markersize=8){

    histo->SetLineColor(color);
    histo->SetLineStyle(linestyle);
    histo->SetLineWidth(linewidth);
    histo->SetMarkerColor(color);
    histo->SetMarkerStyle(markerstyle);
    histo->SetMarkerSize(markersize);

  }

  void FillWithDimColor(TH1* h, bool usealpha=false, float dim=0.8)
  {
    if ( usealpha ){
      h->SetFillColorAlpha(h->GetLineColor(),dim);
      return;
    }
    TColor *color = gROOT->GetColor(h->GetLineColor());
    float R,G,B,hR,hG,hB,hHue,hSat,hVal;
    color->GetRGB(hR,hG,hB);
    color->RGB2HSV(hR,hG,hB,hHue,hSat,hVal);
    color->HSV2RGB(hHue,dim*hSat,hVal,R,G,B);
    h->SetFillColor(color->GetColor(R,G,B));
  }

  // Legends and Texts
  //--------------------------------------------------
  void DrawComponentsLegend(TH1* hnue, TH1* hnumu, TH1* hnc, TH1* hcos, TH1* hother){
  // TLegend *l = new TLegend(0.60, 0.65, 0.85, 0.85, NULL,"brNDC");
    TLegend *l = new TLegend(0.70, 0.70, 0.85, 0.85, NULL,"brNDC");
    l->SetFillStyle(0);
    l->SetTextSize(0.035);
    // l->SetHeader(pot_tag.c_str());
    // l->SetHeader("Integral");
    // l->AddEntry(hnue,   Form("#nu_{e} CC:   %.2f", hnue->Integral()),   "l");
    // l->AddEntry(hnumu,  Form("#nu_{#mu} CC: %.2f", hnumu->Integral()),  "l");
    // l->AddEntry(hnc,    Form("NC:           %.2f", hnc->Integral()),    "l");
    // l->AddEntry(hcos,   Form("Cosmics:      %.2f", hcos->Integral()),   "l");
    // l->AddEntry(hother, Form("Other bkg:    %.2f", hother->Integral()), "l");
    l->AddEntry(hnue,   "#nu_{e} CC",   "l");
    l->AddEntry(hnumu,  "#nu_{#mu} CC", "l");
    l->AddEntry(hnc,    "NC",           "l");
    l->AddEntry(hcos,   "Cosmics",      "l");
    l->AddEntry(hother, "Other bkg",    "l");
    l->Draw("");
  }

  void DrawSigBkgLegend(TH1* h1, char *name1, TH1* h2, char *name2){

    TLegend *leg = new TLegend(.60,.60,.8,.8);
    leg->AddEntry(h1, name1,"l");
    leg->AddEntry(h2, name2,"l");
    leg->SetBorderSize(0); //no border for legend
    leg->SetFillColor(0);  //fill colour is white
    leg->SetFillStyle(0);  //fill colour is white
    leg->SetTextSize(0.04);
    leg->Draw();

  }


  void DrawSigBkgIntLegend(TH1* h1, char *name1, double iSig, TH1* h2, char *name2, double iBkg){

    //double iRatio = iBac/iSig;

    TLegend *leg = new TLegend(.60,.60,.8,.8);
    leg->AddEntry(h1, name1,"l");
    leg->AddEntry(h1, TString::Format("%.2f",iSig),"");
    leg->AddEntry(h2, name2,"l");
    leg->AddEntry(h2, TString::Format("%.2f",iBkg),"");
    //leg->AddEntry(h2, TString::Format("Ratio=%.2f",iRatio),"");
    leg->SetBorderSize(0); //no border for legend                                                                                                                                                                                                                             
    leg->SetFillColor(0);  //fill colour is white                                                                                                                                                                                                                             
    leg->SetFillStyle(0);  //fill colour is white                                                                                                                                                                                                                             
    leg->SetTextSize(0.04);
    leg->Draw();

  }

  //--------------------------------------------------
  void DrawIntEffPurLegend(TH1D* g1, char *name1, TH1D* g2, char *name2){
    
    float int1 = g1->Integral();
    float int2 = g2->Integral();
    TLegend *leg = new TLegend(.65,.2,.85,.3);
    leg->AddEntry(g1, name1,"l");
    leg->AddEntry(g2, name2,"l");
    // leg->AddEntry(g1, Form("%s: %2.f", name1, int1),"l");
    // leg->AddEntry(g2, Form("%s: %2.f", name2, int2),"l");
    leg->SetBorderSize(0); //no border for legend
    leg->SetFillColor(0);  //fill colour is white
    leg->SetFillStyle(0);  //fill colour is white
    leg->SetTextSize(0.03); // this goes in the split canvas so use smaller text
    leg->Draw();

  }

  //--------------------------------------------------
  void DrawPurLegend(TH1D* g1, char *name1){

    TLegend *leg = new TLegend(.65,.25,.85,.35);
    leg->AddEntry(g1, name1,"l");
    leg->SetBorderSize(0); //no border for legend
    leg->SetFillColor(0);  //fill colour is white
    leg->SetFillStyle(0);  //fill colour is white
    leg->SetTextSize(0.03); // this goes in the split canvas so use smaller text
    leg->Draw();

  }

  //--------------------------------------------------
  void DrawEffPurLegend(TH1D* g1, char *name1, TH1D* g2, char *name2){

    TLegend *leg = new TLegend(.65,.2,.85,.3);
    leg->AddEntry(g1, name1,"l");
    leg->AddEntry(g2, name2,"l");
    leg->SetBorderSize(0); //no border for legend
    leg->SetFillColor(0);  //fill colour is white
    leg->SetFillStyle(0);  //fill colour is white
    leg->SetTextSize(0.03); // this goes in the split canvas so use smaller text
    leg->Draw();

  }

  void DrawSigBkgIntText(TH1* hsig, TH1* hbkg, float textsize){

    float isig = hsig->Integral();
    float ibkg = hbkg->Integral();
    float psig = 100. * isig / (isig + ibkg);
    float pbkg = 100. * ibkg / (isig + ibkg);

    TPaveText *pText1 = new TPaveText(0.15, 0.84, 0.30, 0.89, "brNDC");
    TText *text1 = (pText1->AddText("6.6 #times 10^{20} POT"));
    text1->SetTextSize(textsize);
    pText1->SetBorderSize(0);
    pText1->SetFillStyle(0);
    pText1->Draw();
    TPaveText *pText2 = new TPaveText(0.15, 0.72, 0.30, 0.80, "brNDC");
    TText *text2 = pText2->AddText(Form("Sig: %2.f = %2.f %%", isig, psig));
    text2->SetTextAlign(11);
    text2->SetTextSize(textsize);
    TText *text3 = pText2->AddText(Form("Bkg: %2.f = %2.f %%", ibkg, pbkg));
    text3->SetTextAlign(11);
    text3->SetTextSize(textsize);
    pText2->SetBorderSize(0);
    pText2->SetFillStyle(0);
    pText2->Draw();
  }

  // Efficiency and purity graphs 
  //-----------------------------------------------------------------

  TH1D* EffOrPurHistogram(TH1* hSelSignal, TH1* hSelBack, TH1* hSignal, bool geteff) {

    //
    // Make a ROC TGraph for the given signal and background histos
    //
    TH1D* hTotal = (TH1D*)hSelSignal->Clone();
    hTotal->Add(hSelBack);

    TH1D* hPurity = (TH1D*)hSelSignal->Clone();
    hPurity->Divide(hTotal);

    TH1D* hEfficiency = (TH1D*)hSelSignal->Clone();
    hEfficiency->Divide(hSignal);

    if ( geteff )
      return hEfficiency;
    else
      return hPurity;

  }

  TGraph* EffOrPurGraph(TH1* hSelSignal, TH1* hSelBack, TH1* hSignal, bool geteff) {

    const int NBins = hSignal->GetNbinsX();

    TString xTitle = hSignal->GetXaxis()->GetTitle();
    TString yTitle = hSignal->GetYaxis()->GetTitle();

    double eff[NBins], pur[NBins], val[NBins];
    double sb[NBins], ssb[NBins];

    // loop over bins to calculate eff and pur as functions of the bin
    for(unsigned int i = 1; i <= (unsigned int)NBins; ++i) {
      double allsig = hSignal   ->GetBinContent(i);
      double selsig = hSelSignal->GetBinContent(i);
      double selbac = hSelBack  ->GetBinContent(i);

      val[i-1] = hSignal->GetBinCenter(i);
      double EFF = selsig;
      double PUR = selsig;

      if ( (selsig + selbac) > 0 ){
        if ( allsig > 0 ) EFF = EFF/allsig;
        PUR = PUR / (selsig + selbac);
      }

      eff[i-1] = EFF;
      pur[i-1] = PUR;
      if ( PUR > 1 || EFF > 1 ){
        std::cout << "\n\n\n" << " >>> EFF " << EFF << "\t PUR " << PUR << std::endl;
      }
    }

    TString n = hSignal->GetName();
    TGraph *graphPur= new TGraph(NBins,val,pur);
    graphPur->SetName("SelPur_"+n);
    graphPur->GetXaxis()->SetTitle(xTitle);
    graphPur->GetYaxis()->SetTitle("Pur.");
    TGraph *graphEff= new TGraph(NBins,val,eff);
    graphEff->SetName("SelEff_"+n);
    graphEff->GetXaxis()->SetTitle(xTitle);
    graphEff->GetYaxis()->SetTitle("Eff.");

    if ( geteff )
      return graphEff;
    else
      return graphPur;

  }

}
