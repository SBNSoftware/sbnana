#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Ratio.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/UniverseOracle.h"
using namespace ana;

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/TruthCuts.h"

#include "TCanvas.h"
#include "TDecompSVD.h"
#include "TDecompLU.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TPad.h"
#include "TRandom3.h"

template<class T> T sqr(T x){return x*x;}
template<class T> T cube(T x){return x*x*x;}

const double pot = 6e20;

const std::vector<std::string> systs = {
  "DISAth",
  "DISBth",
  "DISCv1u",
  "DISCv2u",
  "IntraNukeNabs",
  "IntraNukeNcex",
  "IntraNukeNinel",
  "IntraNukeNmfp",
  "IntraNukeNpi",
  "IntraNukePIabs",
  "IntraNukePIcex",
  "IntraNukePIinel",
  "IntraNukePImfp",
  "IntraNukePIpi",
  "NC",
  "NonResRvbarp1pi",
  //  "NonResRvbarp1pi", // need to do something with the Alt variants...
  "NonResRvbarp2pi",
  //  "NonResRvbarp2pi",
  "NonResRvp1pi",
  //  "NonResRvp1pi",
  "NonResRvp2pi",
  //  "NonResRvp2pi",
  "ResDecayGamma",
  "CCResAxial",
  "CCResVector",
  //  "CohMA", // return lots of NaNs...
  //  "CohR0",
  "NCELaxial",
  "NCELeta",
  "NCResAxial",
  "NCResVector",
  "QEMA",
};

std::vector<double> local_linear(const std::vector<float>& xs,
                                 const std::vector<double>& ys,
                                 std::vector<double>* grads = 0)
{
  assert(xs.size() == ys.size());

  std::vector<double> ret;

  const double window = 1; // in sigmas

  for(unsigned int i = 0; i < xs.size(); ++i){
    const double x0 = xs[i];

    // Linear fit variables
    double Swx  = 0;
    double Swy  = 0;
    double Swxy = 0;
    double Swy2 = 0;
    double Swx2 = 0;
    double Sw   = 0;

    for(unsigned int j = 0; j < xs.size(); ++j){
      const double x = xs[j];
      const double y = ys[j];
      const double dx = fabs(x-x0)/window;
      if(dx > 1) continue;
      const double w = cube(1-cube(dx));

      Sw   += w;
      Swx  += w*x;
      Swy  += w*y;
      Swx2 += w*x*x;
      Swxy += w*x*y;
      Swy2 += w*y*y;
    } // end for j

    const double d = Sw*Swx2 - Swx*Swx;
    if(d == 0){
      //      std::cout << "D IS ZERO!!!" << std::endl;
      // happens when there is a big gap. Just set the fit to the point itself
      ret.push_back(ys[i]);
      if(grads) grads->push_back(0);
      continue;
    }
    const double m = (Sw*Swxy  - Swx*Swy)/d;
    const double c = (Swy*Swx2 - Swx*Swxy)/d;

    ret.push_back(m*x0+c);
    if(grads) grads->push_back(m);
  } // end for i

  return ret;
}

double calc_mse(const std::vector<std::vector<double>>& preds,
                const std::vector<double>& ys)
{
  const unsigned int Npt = ys.size();
  const unsigned int Nvar = preds.size();

  std::vector<double> p(Npt, 0);
  for(unsigned int jvar = 0; jvar < Nvar; ++jvar){
    for(unsigned int ipt = 0; ipt < Npt; ++ipt) p[ipt] += preds[jvar][ipt];
  }
  double mse = 0;
  for(unsigned int ipt = 0; ipt < Npt; ++ipt) mse += sqr(p[ipt]-ys[ipt]);
  mse /= Npt;
  return mse;
}

// Updates elements of xs[ivar]
std::vector<double> local_linear_update_basis(int ivar,
                                              std::vector<std::vector<float>>& xs,
                                              const std::vector<double>& ys)
{
  const unsigned int Npt = ys.size();
  const unsigned int Nvar = xs.size();

  // try and learn a better beta
  for(int pass = 0; pass < 3; ++pass){ // todo some kind of mse check
    std::vector<double> grads;
    std::vector<double> preds = local_linear(xs[ivar], ys, &grads);
    if(pass == 2/*99*/) return preds;

    // This is taken from
    // https://en.wikipedia.org/wiki/Projection_pursuit_regression#Model_estimation
    TMatrixD W(Npt, Npt);
    for(unsigned int i = 0; i < Npt; ++i) W(i, i) = sqr(grads[i]);
    TVectorD b(Npt);
    TMatrixD X(Npt, Nvar);
    TMatrixD XT(Nvar, Npt);
    for(unsigned int i = 0; i < Npt; ++i){
      b[i] = xs[ivar][i] + (ys[i] - preds[i])/grads[i];
      for(unsigned int j = 0; j < Nvar; ++j){
        X(i, j) = xs[j][i];
        XT(j, i) = xs[j][i];
      }
    }

    bool ok;
    TVectorD beta = TDecompLU(XT*W*X).Solve(XT*W*b, ok);
    if(!ok) return preds;

    // Normalize beta vector
    double norm = 0;
    for(unsigned int j = 0; j < Nvar; ++j) norm += sqr(beta[j]);
    if(isnan(norm) || isinf(norm)) return preds;
    norm = sqrt(norm);
    for(unsigned int j = 0; j < Nvar; ++j) beta[j] /= norm;

    // Re-project the relevant x
    for(unsigned int i = 0; i < Npt; ++i){
      double x = 0;
      for(unsigned int j = 0; j < Nvar; ++j){
        x += beta[j]*xs[j][i];
      }
      xs[ivar][i] = x;
    }
  }

  abort();
}

void plot_residuals(const std::vector<std::vector<double>>& preds,
                    const std::vector<double>& ys)
{
  const unsigned int Npt = ys.size();
  const unsigned int Nvar = preds.size();

  std::vector<double> p(Npt, 0);
  for(unsigned int jvar = 0; jvar < Nvar; ++jvar){
    for(unsigned int ipt = 0; ipt < Npt; ++ipt) p[ipt] += preds[jvar][ipt];
  }
  TGraph* g = new TGraph;
  for(unsigned int ipt = 0; ipt < Npt; ++ipt)
    g->SetPoint(ipt, p[ipt], ys[ipt]-p[ipt]);

  g->SetMarkerStyle(kFullDotMedium);
  g->Draw("ap");
}

void plot_preds(const std::vector<std::vector<float>>& xs,
                const std::vector<double>& ys,
                const std::vector<std::vector<double>>& preds)
{
  const unsigned int Npt = ys.size();
  const unsigned int Nvar = xs.size();

  TCanvas* c = new TCanvas;
  unsigned int Nx = 0, Ny = 0;
  while(Nx*Ny < Nvar) if(Nx > Ny) ++Ny; else ++Nx;
  c->Divide(Nx, Ny);

  for(unsigned int ivar = 0; ivar < Nvar; ++ivar){
    TGraph* gdat = new TGraph;
    TGraph* gpred = new TGraph;

    for(unsigned int ipt = 0; ipt < Npt; ++ipt){
      double p = 0;
      for(unsigned int jvar = 0; jvar < Nvar; ++jvar){
        if(jvar == ivar) continue;
        p += preds[jvar][ipt];
      }

      gdat->SetPoint(gdat->GetN(), xs[ivar][ipt], ys[ipt]-p);
      gpred->SetPoint(gpred->GetN(), xs[ivar][ipt], preds[ivar][ipt]);
    }

    c->cd(ivar+1);
    gdat->SetMarkerStyle(kFullDotMedium);
    gdat->Draw("ap");
    gpred->Sort();
    gpred->SetLineColor(kRed);
    gpred->Draw("l same");
  }

  c->cd(0);
}

std::vector<std::vector<double>> additive_model(/*const*/ std::vector<std::vector<float>>& xs,
                                                const std::vector<double>& ys)
{
  const unsigned int Npt = ys.size();
  for(const std::vector<float>& x: xs) assert(x.size() == Npt);
  const unsigned int Nvar = xs.size();

  std::cout << "Solving model for " << Npt << " universes described by " << Nvar << " vars" << std::endl;

  std::vector<std::vector<double>> preds(Nvar, std::vector<double>(Npt, 0));

  double old_mse = calc_mse(preds, ys);
  std::cout << "MSE " << old_mse << " (" << sqrt(old_mse) << ")" << std::endl;

  plot_preds(xs, ys, preds);
  gPad->Print("preds_orig.pdf");

  for(int pass = 0; pass < 100; ++pass){
    for(unsigned int ivar = 0; ivar < Nvar; ++ivar){
      // Residual
      std::vector<double> dy = ys;
      for(unsigned int jvar = 0; jvar < Nvar; ++jvar){
        if(jvar == ivar) continue; // Not including this variable
        for(unsigned int ipt = 0; ipt < Npt; ++ipt) dy[ipt] -= preds[jvar][ipt];
      } // end for jvar

      preds[ivar] = local_linear(xs[ivar], dy);
    } // end for ivar

    plot_preds(xs, ys, preds);
    gPad->Print(TString::Format("preds_%d.pdf", pass).Data());

    const double mse = calc_mse(preds, ys);
    std::cout << pass << ": MSE " << mse << " (" << sqrt(mse) << ")" << std::endl;

    if(mse >= old_mse) break; // convergence
    old_mse = mse;
  } // end for pass

  for(int pass = 0; pass < 100; ++pass){
    for(unsigned int ivar = 0; ivar < Nvar; ++ivar){
      // Residual
      std::vector<double> dy = ys;
      for(unsigned int jvar = 0; jvar < Nvar; ++jvar){
        if(jvar == ivar) continue; // Not including this variable
        for(unsigned int ipt = 0; ipt < Npt; ++ipt) dy[ipt] -= preds[jvar][ipt];
      } // end for jvar

      preds[ivar] = local_linear_update_basis(ivar, xs, dy);
    } // end for ivar

    plot_preds(xs, ys, preds);
    gPad->Print(TString::Format("preds_%d.pdf", pass).Data());

    const double mse = calc_mse(preds, ys);
    std::cout << pass << ": MSE " << mse << " (" << sqrt(mse) << ")" << std::endl;

    if(mse >= old_mse) break; // convergence
    old_mse = mse;
  } // end for pass

  plot_residuals(preds, ys);
  gPad->Print("residuals.pdf");

  //  plot_residuals(preds, ys);
  //  gPad->Print("residuals.pdf");

  return preds;
}

void test_additive(bool reload = false)
{
  if(reload || TFile("additive_state.root").IsZombie()){
    SpectrumLoader loader("workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_intrnue_spill_gsimple-configf-v1_tpc_flat_caf_sbnd");

    HistAxis axis("True Energy (GeV)", Binning::Simple(20, 0, 5), SIMPLEVAR(truth.E));

    Spectrum snom(loader, axis, kIsCC);

    std::vector<Spectrum> multiverse;
    multiverse.reserve(100);
    for(int i = 0; i < 100; ++i) multiverse.emplace_back(loader, axis, kIsCC, kNoShift, GetUniverseWeight(systs, i));

    loader.Go();

    const UniverseOracle& uo = UniverseOracle::Instance();
    std::vector<std::vector<float>> xs;
    for(const std::string& s: systs) xs.push_back(uo.ShiftsForSyst(s));

    TFile fout("additive_state.root", "RECREATE");
    snom.SaveTo(fout.mkdir("nominal"));
    for(int i = 0; i < 100; ++i) multiverse[i].SaveTo(fout.mkdir(TString::Format("multiverse_%d", i).Data()));
    for(unsigned int i = 0; i < systs.size(); ++i){
      TVectorD v(100);
      for(int j = 0; j < 100; ++j) v[j] = xs[i][j];
      v.Write(TString::Format("xs_%d", i).Data());
    }
  }

  TFile* fin = new TFile("additive_state.root");
  Spectrum snom = *LoadFrom<Spectrum>(fin->GetDirectory("nominal"));
  std::vector<Spectrum> multiverse;
  for(int i = 0; i < 100; ++i) multiverse.push_back(*LoadFrom<Spectrum>(fin->GetDirectory(TString::Format("multiverse_%d", i).Data())));

  std::vector<std::vector<float>> xs;
  for(unsigned int i = 0; i < systs.size(); ++i){
    TVectorD* v = (TVectorD*)fin->Get(TString::Format("xs_%d", i).Data());
    xs.emplace_back();
    for(int j = 0; j < 100; ++j) xs.back().push_back((*v)[j]);
  }

  std::vector<double> ys;

  for(int i = 0; i < 100; ++i){
    // TODO leaks histogram
    const double y = (multiverse[i] / snom).ToTH1()->GetBinContent(5);
    if(isnan(y) || isinf(y) || isnan(log(y)) || isinf(log(y))){
      std::cout << "bad y = " << y << std::endl;
      continue;
    }

    ys.push_back(log(y));
  }

  // Example of switching to a different basis - gives worse results if you
  // just use a random one...
  std::vector<TVectorD> basis(systs.size(), TVectorD(systs.size()));

  // Random matrix with unit vector columns
  for(unsigned int i = 0; i < systs.size(); ++i){
    double norm = 0;
    for(unsigned int j = 0; j < systs.size(); ++j){
      basis[i][j] = gRandom->Gaus();
      norm += sqr(basis[i][j]);
    }
    norm = sqrt(norm);
    for(unsigned int j = 0; j < systs.size(); ++j){
      basis[i][j] /= norm;
    }
  }


  std::vector<std::vector<float>> newxs;
  newxs.resize(xs.size());
  for(auto& it: newxs) it.resize(xs[0].size());

  for(int i = 0; i < 100; ++i){ // for each universe
    // Work out the new x-values by dotting with the basis
    for(unsigned int j = 0; j < systs.size(); ++j){
      newxs[j][i] = 0;
      for(unsigned int k = 0; k < systs.size(); ++k){
        newxs[j][i] += xs[k][i] * basis[j][k];
      }
    }
  }

  // Use the new basis instead
  //  xs = newxs;
  // Use both natural and new basis
  //  xs.insert(xs.end(), newxs.begin(), newxs.end());

  /*
  TH1* h = snom.ToTH1(pot);
  h->Draw("hist");
  h->GetYaxis()->SetRangeUser(0, 1.3*h->GetMaximum());
  for(const Spectrum& s: multiverse) s.ToTH1(pot, kGray)->Draw("hist same");
  h->Draw("hist same");
  gPad->Print("multiverse.pdf");
  */

  additive_model(xs, ys);
}
