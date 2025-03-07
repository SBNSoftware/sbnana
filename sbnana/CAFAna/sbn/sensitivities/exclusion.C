#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Vars/FitVarsSterileApprox.h"
#include "sbnana/CAFAna/Prediction/PredictionInterp.h"
#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"
#include "sbnana/CAFAna/Experiment/MultiExperimentSBN.h"
#include "sbnana/CAFAna/Experiment/CountingExperiment.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/Analysis/Surface.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
using namespace ana;

#include "OscLib/IOscCalc.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraph.h"

#include <fstream>
#include <vector>

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void exclusion(const std::string anatype = numuStr)
{
  const char* name_in;
  const char* name_out;
  if (anatype == numuStr) {
    name_in = "surfaces_numu.root";
    name_out = "exclusion_numu_quick_count.pdf";
  }
  else if (anatype == nueStr) {
    name_in = "surfaces_nue.root";
    name_out = "exclusion_nue.pdf";
  }
  else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }

  TFile fin(name_in);  

  std::vector<std::pair<std::string, std::string>> selections{
    {"", "Fake Reco"},
    //{"_np", "1#muNp"},
    //{"_1p", "1#mu1p"}
  };

  std::vector<std::vector<std::tuple<std::string, std::string, Color_t>>> plots{
    {
      {"xsec only", "_xsec", kBlack}
    },
    {
      {"flux only", "_flux", kBlack}
    },
    /*{
      {"2% normalization only", "_norm", kBlack}
    },
    {
      {"10% normalization only", "_big_norm", kBlack}
    },
    {
      {"energy scale only", "_energy", kBlack}
    },*/
    {
      {"xsec + flux", "_xsec_flux", kBlack},
      //{"xsec + flux + 2% norm", "_xsec_flux_norm", kBlue},
      //{"xsec + flux + 2% norm + energy", "_xsec_flux_norm_energy", kRed}
    },
    //{
    //  {"xsec + flux", "_xsec_flux", kBlack},
    //  {"xsec + flux + 10% norm", "_xsec_flux_big_norm", kBlue},
    //  {"xsec + flux + 10% norm + energy", "_xsec_flux_big_norm_energy", kRed}
    //}
  }; 

  TCanvas c;
  c.Print((name_out + "["s).c_str());

  for(const auto& [sel_suffix, sel_name]: selections) {
    Surface& surf_nom = *ana::LoadFrom<Surface>(fin.GetDirectory(("exclusion/nom_fd"+sel_suffix).c_str())).release();

    TH2* crit5sig = Gaussian5Sigma1D1Sided(surf_nom);
    TH2* crit3sig = Gaussian3Sigma1D1Sided(surf_nom);
    TH2* crit90 = Gaussian90Percent1D1Sided(surf_nom);
    TH2* crit95 = Gaussian95Percent1D1Sided(surf_nom);
    TH2* crit99 = Gaussian99Percent1D1Sided(surf_nom);

    std::vector<std::tuple<std::string, TH2*>> sigs{
      {"5#sigma", crit5sig},
      {"3#sigma", crit3sig},
      {"99%", crit99},
      {"95%", crit95},
      {"90%", crit90}
    };
    
    for(const auto& [sig_name, sig]: sigs) {
      for(const auto& curves: plots) {
        surf_nom.SetTitle((sig_name+" Exclusion - " + sel_name).c_str());
        surf_nom.DrawContour(sig, 7, kBlack);

        TLegend * lgdis = new TLegend(0.15,0.15,0.45,0.45);
        lgdis->SetFillColor(0);
        lgdis->SetBorderSize(0);
        TH1* dummy = new TH1F("", "", 1, 0, 1);

        for(const auto& [leg, ext, col]: curves) {
          Surface& surf_syst = *ana::LoadFrom<Surface>(fin.GetDirectory(("exclusion/fd"+ext+sel_suffix).c_str())).release(); 
          surf_syst.DrawContour(sig, kSolid, col);
          dummy->SetLineColor(col);
          lgdis->AddEntry(dummy->Clone(), leg.c_str());
        }
        dummy->SetLineStyle(7);
        dummy->SetLineColor(kBlack);
        lgdis->AddEntry(dummy->Clone(), "Stats Only");
        lgdis->Draw("same");

        c.Print(name_out);
        c.Clear();
      }
    }
  }

  c.Print((name_out + "]"s).c_str());
    
  //proposal_90pctCL->SetLineColor(kGreen);
  //proposal_90pctCL->Draw("l");
  //minosp_90pctCL->SetLineColor(kCyan);
  //minosp_90pctCL->Draw("l");
  //minisciboone_90pctCL->SetLineColor(kOrange);
  //minisciboone_90pctCL->Draw("l");

  //lgdis->AddEntry(proposal_90pctCL, "Proposal 90%");
  //lgdis->AddEntry(minosp_90pctCL, "Minos/Minos+ 90%");
  //lgdis->AddEntry(minisciboone_90pctCL, "MiniBoone/SciBoone 90%");

//  TFile fout("exclusion_graphs.root", "RECREATE");
//  
//  std::vector<TGraph*> graph_nd_fd_90 = surf_nom_nd_fd.GetGraphs(crit90);
//  std::vector<TGraph*> graph_nd_fd_95 = surf_nom_nd_fd.GetGraphs(crit95);
//  std::vector<TGraph*> graph_nd_fd_99 = surf_nom_nd_fd.GetGraphs(crit99);
//  std::vector<TGraph*> graph_nd_fd_3s = surf_nom_nd_fd.GetGraphs(crit3sig);
//  std::vector<TGraph*> graph_nd_fd_5s = surf_nom_nd_fd.GetGraphs(crit5sig);
//  
//  std::vector<Surface> surfaces {surf_nom, surf_nom_nd, surf_nom_ub, surf_nom_fd, surf_nom_nd_fd, surf_syst_all, surf_syst_nd, surf_syst_ub, surf_syst_fd, surf_syst_nd_fd};
//
//  std::vector<TH2*> con_lev {crit90, crit95, crit99, crit3sig, crit5sig};
//
//  std::vector<std::string> s_name {"all_stat", "nd_stat", "ub_stat", "fd_stat", "nd_fd_stat", "all_syst", "nd_syst", "ub_syst", "fd_syst", "nd_fd_syst"};
//
//  std::vector<std::string> c_name {"_90pct", "_95pct", "_99pct", "_3sig", "_5sig"};
//  
//  std::cout<<"Before Loop"<<std::endl;
// 
//  ofstream txtfile;
//  txtfile.open("CAFAna_exclusion.txt");
//
//  //std::vector<double> x_v;
//  //std::vector<double> y_v;
//
//  for (int i = 0; i < 10; ++i) {
//    txtfile<<s_name[i]<<"\n"<<"=========="<<"\n";
//    for (int j = 0; j < 5; ++j) {
//      std::vector<TGraph*> graphs = surfaces[i].GetGraphs(con_lev[j]);
//      std::cout<<"Got Graph"<<std::endl;
//      txtfile<<c_name[j]<<"\n"<<"------"<<"\n";
//      for (auto g : graphs) {
//        g->Write((s_name[i]+c_name[j]).c_str());
//        int n = g->GetN();
//        double *x = g->GetX();
//        double *y = g->GetY();
//        for (int k = 0; k < n; ++k) {
//	  //x_v.push_back(x[k]);
//	  //y_v.push_back(y[k]);
//          txtfile<<x[k]<<" "<<y[k]<<"\n";
//	}
//      }
//      txtfile<<"\n";
//    }
//    txtfile<<"\n\n";
//  } 
//
//  txtfile.close();

}
