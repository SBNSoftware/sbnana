#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Vars/FitVarsSterileApprox.h"
#include "sbnana/CAFAna/Experiment/IExperiment.h"
#include "sbnana/CAFAna/Experiment/MultiExperiment.h"
#include "sbnana/CAFAna/Analysis/Surface.h"

#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"

using namespace ana;

void contours_overlay()
{
  // contours_icarus.C
  IExperiment* expt_fd = LoadFromFile<IExperiment>("icarus_numu_disapp_exclusion_pac2021.root", "fd_expt").release();

  // contours_sbnd.C
  IExperiment* expt_nd = LoadFromFile<IExperiment>("sbnd_numu_disapp_exclusion_pac2021.root", "nd_expt").release();

  MultiExperiment expt_joint({expt_fd, expt_nd});

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  const FitAxis kAxTh(&kFitSinSq2ThetaMuMu, 60, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 60, 1e-2, 1e2, true);

  Surface surf_fd(expt_fd, calc, kAxTh, kAxDmSq);
  Surface surf_nd(expt_nd, calc, kAxTh, kAxDmSq);
  Surface surf_joint(&expt_joint, calc, kAxTh, kAxDmSq);

  TH2* crit5sig = Gaussian5Sigma1D1Sided(surf_fd);
  surf_joint.DrawContour(crit5sig, kSolid, kGreen+2);
  surf_fd.DrawContour(crit5sig, kSolid, kRed);
  surf_nd.DrawContour(crit5sig, kSolid, kBlue);

  // slack message from Jacob in #oscillation on 5 May 2021
  TFile fin("CAFAna_numu_disapp_exclusion.root");
  TGraph* g_fd = (TGraph*)fin.Get("fd_stat_5sig");
  g_fd->SetLineColor(kRed);
  g_fd->SetLineWidth(1);
  g_fd->DrawClone("l same");

  TGraph* g_nd = (TGraph*)fin.Get("nd_stat_5sig");
  g_nd->SetLineColor(kBlue);
  g_nd->SetLineWidth(1);
  g_nd->DrawClone("l same");

  TGraph* g_joint = (TGraph*)fin.Get("nd_fd_stat_5sig");
  g_joint->SetLineColor(kGreen+2);
  g_joint->SetLineWidth(1);
  g_joint->DrawClone("l same");

  TLegend* leg = new TLegend(.125, .125, .55, .5);
  leg->SetFillStyle(0);
  leg->AddEntry((TObject*)0, "Stat only, 5#sigma C.L.", "");
  leg->AddEntry((TObject*)0, "6.6#times10^{20} POT", "");
  TH1* dummy = new TH1F("", "", 1, 0, 1);
  dummy->SetLineColor(kRed);
  leg->AddEntry(dummy->Clone(), "Icarus", "l");
  dummy->SetLineColor(kBlue);
  leg->AddEntry(dummy->Clone(), "SBND", "l");
  dummy->SetLineColor(kGreen+2);
  leg->AddEntry(dummy->Clone(), "SBND+Icarus", "l");
  dummy->SetLineColor(kBlack);
  dummy->SetLineWidth(1);
  leg->AddEntry(dummy->Clone(), "Parameterized", "l");
  dummy->SetLineWidth(2);
  leg->AddEntry(dummy->Clone(), "Full reco.", "l");
  leg->Draw();

  gPad->Print("contour_overlay.pdf");
  gPad->Print("contour_overlay.png");

  TH2* crit3sig = Gaussian3Sigma1D1Sided(surf_fd);
  TH2* crit90pc = Gaussian90Percent1D1Sided(surf_fd);
  TH2* crit95pc = Gaussian95Percent1D1Sided(surf_fd);
  TH2* crit99pc = Gaussian99Percent1D1Sided(surf_fd);

  TFile fout("joint_numu_disapp_exclusion_pac2021.root", "RECREATE");
  surf_joint.GetGraphs(crit3sig)[0]->Write("nd_fd_stat_3sig");
  surf_joint.GetGraphs(crit5sig)[0]->Write("nd_fd_stat_5sig");
  surf_joint.GetGraphs(crit90pc)[0]->Write("nd_fd_stat_90pct");
  surf_joint.GetGraphs(crit95pc)[0]->Write("nd_fd_stat_95pct");
  surf_joint.GetGraphs(crit99pc)[0]->Write("nd_fd_stat_99pct");
}
