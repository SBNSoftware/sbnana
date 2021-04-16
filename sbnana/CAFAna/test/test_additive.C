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

// Returns new predictions
TVectorD local_linear(const TVectorD& xs,
                      const TVectorD& ys,
                      TVectorD* invgrads = 0)
{
  assert(xs.GetNrows() == ys.GetNrows());
  const unsigned int Npts = xs.GetNrows();

  if(invgrads) invgrads->ResizeTo(Npts);

  TVectorD ret(Npts);

  const double window = 1; // in sigmas

  for(unsigned int i = 0; i < Npts; ++i){
    const double x0 = xs[i];

    // Linear fit variables
    double Swx  = 0;
    double Swy  = 0;
    double Swxy = 0;
    double Swy2 = 0;
    double Swx2 = 0;
    double Sw   = 0;

    for(unsigned int j = 0; j < Npts; ++j){
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
      //      std::cout << i << std::endl;
      //      std::cout << "D IS ZERO!!!" << std::endl;
      // happens when there is a big gap. Just set the fit to the point itself
      ret[i] = ys[i];
      if(invgrads) (*invgrads)[i] = std::numeric_limits<double>::infinity(); // TODO what should we do?
      continue;
    }

    const double m = (Sw*Swxy  - Swx*Swy)/d;
    const double c = (Swy*Swx2 - Swx*Swxy)/d;

    ret[i] = m*x0+c;
    if(invgrads) (*invgrads)[i] = 1/m;
  } // end for i

  return ret;
}

TVectorD total_prediction(const std::vector<TVectorD>& preds)
{
  TVectorD tot(preds[0].GetNrows());
  for(const TVectorD& pred: preds) tot += pred;
  return tot;
}

double calc_mse(const std::vector<TVectorD>& preds,
                const TVectorD& ys)
{
  const unsigned int Npt = ys.GetNrows();

  return (total_prediction(preds) - ys).Norm2Sqr() / Npt;
}

TVectorD project(const TVectorD& beta, const std::vector<TVectorD>& xs)
{
  const int Nvar = beta.GetNrows();
  const int Npt = xs[0].GetNrows();

  TVectorD proj(Npt);
  for(int i = 0; i < Npt; ++i){
    for(int j = 0; j < Nvar; ++j){
      proj[i] += beta[j]*xs[j][i];
    }
  }

  return proj;
}

std::vector<TVectorD> project(const std::vector<TVectorD>& betas,
                              const std::vector<TVectorD>& xs)
{
  std::vector<TVectorD> ret;
  ret.reserve(betas.size());
  for(const TVectorD& beta: betas) ret.push_back(project(beta, xs));
  return ret;
}

// Updates beta, returns new predictions
TVectorD local_linear_update_basis(TVectorD& beta,
                                   const std::vector<TVectorD>& xs,
                                   const TVectorD& ys)
{
  const unsigned int Npt = ys.GetNrows();
  const unsigned int Nvar = xs.size();

  double old_mse = std::numeric_limits<double>::infinity();

  TVectorD oldpreds(Npt);
  TVectorD oldbeta(beta.GetNrows());

  // try and learn a better beta
  while(true){
    TVectorD bx = project(beta, xs);
    TVectorD invgrads;
    TVectorD preds = local_linear(bx, ys, &invgrads);

    const double mse = (preds-ys).Norm2Sqr() / Npt;
    if(mse > old_mse){
      //      std::cout << mse << " " << old_mse << " bail after " << pass << std::endl;
      beta = oldbeta;
      return oldpreds;
    }
    old_mse = mse;
    oldbeta = beta;
    oldpreds = preds;

    // This is taken from
    // https://en.wikipedia.org/wiki/Projection_pursuit_regression#Model_estimation
    TMatrixD W(Npt, Npt);
    for(unsigned int i = 0; i < Npt; ++i){
      if(!isinf(invgrads[i]) && !isnan(invgrads[i])) W(i, i) = sqr(1/invgrads[i]);
    }
    TVectorD b(Npt);
    TMatrixD X(Npt, Nvar);
    TMatrixD XT(Nvar, Npt);

    for(unsigned int i = 0; i < Npt; ++i){
      // TODO could almost write in one expression, but we don't have element-wise product
      b[i] = bx[i];
      if(!isinf(invgrads[i]) && !isnan(invgrads[i])) b[i] += (ys[i] - preds[i])*invgrads[i];
      // TODO figure out what to do with zero gradient case!
      //      if(grads[i] != 0) b(i, 0) += (ys[i] - preds[i])/grads[i]; else b(i, 0) += 1e8;
      //      b(i, 0) += (ys[i] - preds[i])/grads[i];
      for(unsigned int j = 0; j < Nvar; ++j){
        X(i, j) = xs[j][i];
        XT(j, i) = xs[j][i];
      }
    }

    bool ok;
    TVectorD newbeta = TDecompLU(XT*W*X).Solve(XT*W*b, ok);
    if(!ok) return preds;

    // Normalize beta vector
    double norm = newbeta.Norm2Sqr();
    if(isnan(norm) || isinf(norm)) return preds;/*{
      W.Print();
      b.Print();
      X.Print();
      XT.Print();
      XWX.Print();
      beta.Print();
      abort();
      return preds; // bail out
      }*/
    norm = sqrt(norm);
    newbeta *= 1/norm;
    beta = newbeta;
  }

  abort(); // unreached
}

void plot_residuals(const std::vector<TVectorD>& preds,
                    const TVectorD& ys)
{
  const unsigned int Npt = ys.GetNrows();
  const unsigned int Nvar = preds.size();

  TVectorD p = total_prediction(preds);
  TGraph* g = new TGraph;
  for(unsigned int ipt = 0; ipt < Npt; ++ipt)
    g->SetPoint(ipt, p[ipt], ys[ipt]-p[ipt]);

  g->SetMarkerStyle(kFullDotMedium);
  g->Draw("ap");
}

void plot_preds(const std::vector<TVectorD>& xs,
                const TVectorD& ys,
                const std::vector<TVectorD>& preds)
{
  const unsigned int Npt = ys.GetNrows();
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

std::vector<TVectorD> additive_model(/*const*/ std::vector<TVectorD>& xs,
                                     const TVectorD& ys)
{
  const unsigned int Npt = ys.GetNrows();
  for(const TVectorD& x: xs) assert(x.GetNrows() == int(Npt));
  const unsigned int Nvar = xs.size();

  std::cout << "Solving model for " << Npt << " universes described by " << Nvar << " vars" << std::endl;

  std::vector<TVectorD> preds(Nvar, TVectorD(Npt));

  double old_mse = calc_mse(preds, ys);
  std::cout << "MSE " << old_mse << " (" << sqrt(old_mse) << ")" << std::endl;

  plot_preds(xs, ys, preds);

  gPad->Print("preds_anim.pdf[");

  gPad->Print("preds_orig.pdf");
  gPad->Print("preds_anim.pdf");
  /*
  for(int pass = 0; pass < 100; ++pass){
    for(unsigned int ivar = 0; ivar < Nvar; ++ivar){
      // Residual
      const TVectorD dy = ys - (total_prediction(preds) - preds[ivar]);
      preds[ivar] = local_linear(xs[ivar], dy);
    } // end for ivar

    plot_preds(xs, ys, preds);
    gPad->Print(TString::Format("preds_%d.pdf", pass).Data());
    gPad->Print("preds_anim.pdf");

    const double mse = calc_mse(preds, ys);
    std::cout << pass << ": MSE " << mse << " (" << sqrt(mse) << ")" << std::endl;

    if(mse >= old_mse) break; // convergence
    old_mse = mse;
  } // end for pass
  */

  std::vector<TVectorD> betas(Nvar, TVectorD(Nvar));
  // Start with the same basis as the regular variables
  for(unsigned int i = 0; i < Nvar; ++i) betas[i][i] = 1;

  /* // random initialization
  for(unsigned int i = 0; i < Nvar; ++i){
    for(unsigned int j = 0; j < Nvar; ++j){
      betas[i][j] = gRandom->Gaus();
    }
    betas[i] *= 1/sqrt(betas[i].Norm2Sqr()); // unit vector
  }
  */

  int pass = 0;
  while(true){
    for(unsigned int ivar = 0; ivar < Nvar; ++ivar){
      // Residual
      const TVectorD dy = ys - (total_prediction(preds) - preds[ivar]);
      preds[ivar] = local_linear_update_basis(betas[ivar], xs, dy);
    } // end for ivar

    plot_preds(project(betas, xs), ys, preds);
    gPad->Print(TString::Format("preds_%d.pdf", pass).Data());
    gPad->Print("preds_anim.pdf");

    const double mse = calc_mse(preds, ys);
    std::cout << pass << ": MSE " << mse << " (" << sqrt(mse) << ")" << std::endl;

    if(mse >= old_mse) break; // convergence
    old_mse = mse;

    ++pass;
  } // end for pass

  gPad->Print("preds_anim.pdf]");

  //  for(const TVectorD& beta: betas) beta.Print();

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

  std::vector<TVectorD> xs;
  for(unsigned int i = 0; i < systs.size(); ++i){
    TVectorD* v = (TVectorD*)fin->Get(TString::Format("xs_%d", i).Data());
    xs.emplace_back(100);
    for(int j = 0; j < 100; ++j) xs.back()[j] = (*v)[j];
  }

  TVectorD ys(100);

  for(int i = 0; i < 100; ++i){
    // TODO leaks histogram
    const double y = (multiverse[i] / snom).ToTH1()->GetBinContent(5);
    if(isnan(y) || isinf(y) || isnan(log(y)) || isinf(log(y))){
      std::cout << "bad y = " << y << std::endl;
      abort();
    }

    /*
    // Out of range in first pass
    if(i == 1 || i == 20 || i == 54 || i == 74 || i == 75 || i == 83 || i == 86 || i == 90){
      ys[i] = 0;
      for(auto& x: xs) x[i] = 0;
      continue;
    }
    */

    /*
    if(log(y) > 0.4 | log(y) < -0.4){
      // Very specific...
      std::cout << "bad y = " << y << std::endl;
      for(auto& x: xs) x.erase(x.begin()+i);
      continue;
    }
    */

    ys[i] = log(y);
  }

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
