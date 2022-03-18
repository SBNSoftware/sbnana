#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Ratio.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Systs/UniverseOracle.h"
using namespace ana;

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/TruthCuts.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TPad.h"
#include "TRandom3.h"

#include <Eigen/Dense>

template<class T> T sqr(T x){return x*x;}
template<class T> T cube(T x){return x*x*x;}

const double pot = 6e20;

const int Ndim = 2; // dimensions in the model

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
Eigen::VectorXd local_linear(const Eigen::VectorXd& xs,
                             const Eigen::VectorXd& ys,
                             const Eigen::VectorXd& xs_test,
                             Eigen::VectorXd* grads = 0)
{
  assert(xs.size() == ys.size());
  const unsigned int Npts = xs.size();
  const unsigned int Ntest = xs_test.size();

  if(grads){
    grads->resize(Ntest);
    grads->setZero();
  }

  Eigen::VectorXd ret = Eigen::VectorXd::Zero(Ntest);

  const double window = 2;//1; // in sigmas

  for(unsigned int i = 0; i < Ntest; ++i){
    const double x0 = xs_test[i];

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
      ret[i] = ys[i]; // TODO not right for test case
      if(grads) (*grads)[i] = 0; // TODO what should we do?
      continue;
    }

    const double m = (Sw*Swxy  - Swx*Swy)/d;
    const double c = (Swy*Swx2 - Swx*Swxy)/d;

    ret[i] = m*x0+c;
    if(grads) (*grads)[i] = m;
  } // end for i

  return ret;
}

Eigen::VectorXd total_prediction(const Eigen::MatrixXd& preds)
{
  // Sum up the contents of each row
  return preds.rowwise().sum();
}

double calc_mse(const Eigen::MatrixXd& preds,
                const Eigen::VectorXd& ys)
{
  const unsigned int Npt = ys.size();

  return (total_prediction(preds) - ys).squaredNorm() / Npt;
}

Eigen::VectorXd projectSingle(const Eigen::VectorXd& beta, const Eigen::MatrixXd& xs)
{
  return xs*beta;
}

Eigen::MatrixXd projectMulti(const Eigen::MatrixXd& betas,
                             const Eigen::MatrixXd& xs)
{
  return xs*betas;
}

// Updates beta, returns new predictions
Eigen::VectorXd local_linear_update_basis(Eigen::VectorXd& beta,
                                          const Eigen::MatrixXd& xs,
                                          const Eigen::VectorXd& ys,
                                          const Eigen::MatrixXd& xs_test,
                                          Eigen::VectorXd& preds_test)
{
  const unsigned int Npt = ys.size();
  const unsigned int Nvar = xs.cols();
  const unsigned int Ntest = xs_test.cols();

  double old_mse = std::numeric_limits<double>::infinity();

  Eigen::VectorXd oldpreds = Eigen::VectorXd::Zero(Npt);
  Eigen::VectorXd oldpreds_test = Eigen::VectorXd::Zero(Ntest);
  Eigen::VectorXd oldbeta = Eigen::VectorXd::Zero(beta.size());

  // try and learn a better beta
  while(true){
    const Eigen::VectorXd bx = projectSingle(beta, xs);

    Eigen::VectorXd grads;
    const Eigen::VectorXd preds = local_linear(bx, ys, bx, &grads);

    const Eigen::VectorXd bx_test = projectSingle(beta, xs_test);
    preds_test = local_linear(bx, ys, bx_test, 0);

    const double mse = (preds-ys).squaredNorm() / Npt;
    if(mse > old_mse){
      //      std::cout << mse << " " << old_mse << " bail after " << pass << std::endl;
      beta = oldbeta;
      preds_test = oldpreds_test;
      return oldpreds;
    }
    old_mse = mse;
    oldbeta = beta;
    oldpreds = preds;
    oldpreds_test = preds_test;

    // This is taken from
    // https://en.wikipedia.org/wiki/Projection_pursuit_regression#Model_estimation
    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(Npt, Npt);
    for(unsigned int i = 0; i < Npt; ++i){
      if(!isinf(grads[i]) && !isnan(grads[i])) W(i, i) = sqr(grads[i]);
    }

    // TODO figure out what to do with zero gradient case!

    // Conversions ensure element-wise operation
    const Eigen::VectorXd b = bx + ((ys-preds).array() / grads.array()).matrix();

    //    const Eigen::VectorXd newbeta = (xs.transpose()*W*xs).colPivHouseholderQr().solve(xs.transpose()*W*b);
    const Eigen::VectorXd newbeta = (xs.transpose()*W*xs).ldlt().solve(xs.transpose()*W*b);
    // TODO detect solution failure and just return existing preds

    // Normalize beta vector
    const double norm = newbeta.squaredNorm();
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

    beta = newbeta/sqrt(norm);
  }

  abort(); // unreached
}

void plot_residuals(const Eigen::MatrixXd& preds,
                    const Eigen::VectorXd& ys)
{
  const unsigned int Npt = ys.size();

  Eigen::VectorXd p = total_prediction(preds);
  TGraph* g = new TGraph;
  for(unsigned int ipt = 0; ipt < Npt; ++ipt)
    g->SetPoint(ipt, p[ipt], ys[ipt]-p[ipt]);

  g->SetMarkerStyle(kFullDotMedium);
  g->Draw("ap");
}

void plot_preds(const Eigen::MatrixXd& xs,
                const Eigen::VectorXd& ys,
                const Eigen::MatrixXd& preds,
                const Eigen::MatrixXd& xs_test,
                const Eigen::VectorXd& ys_test,
                const Eigen::MatrixXd& preds_test)
{
  const unsigned int Npt = ys.size();
  const unsigned int Ntest = ys_test.size();
  const unsigned int Nvar = xs.cols();

  TCanvas* c = new TCanvas;
  unsigned int Nx = 0, Ny = 0;
  while(Nx*Ny < Ndim) if(Nx > Ny) ++Ny; else ++Nx;
  c->Divide(Nx, Ny);

  for(unsigned int idim = 0; idim < Ndim; ++idim){
    TGraph* gdat = new TGraph;
    TGraph* gdat_test = new TGraph;
    TGraph* gpred = new TGraph;

    for(unsigned int ipt = 0; ipt < Npt; ++ipt){
      double p = 0;
      for(unsigned int jdim = 0; jdim < Ndim; ++jdim){
        if(jdim == idim) continue;
        p += preds(ipt, jdim);
      }

      gdat->SetPoint(gdat->GetN(), xs(ipt, idim), ys[ipt]-p);
      gpred->SetPoint(gpred->GetN(), xs(ipt, idim), preds(ipt, idim));
    }

    for(unsigned int itest = 0; itest < Ntest; ++itest){
      double p = 0;
      for(unsigned int jdim = 0; jdim < Ndim; ++jdim){
        if(jdim == idim) continue;
        p += preds_test(itest, jdim);
      }

      gdat_test->SetPoint(gdat_test->GetN(), xs_test(itest, idim), ys_test[itest]-p);
      // Cross-check not using crazy predictions
      gpred->SetPoint(gpred->GetN(), xs_test(itest, idim), preds_test(itest, idim));
    }

    c->cd(idim+1);
    gdat->SetMarkerStyle(kFullDotMedium);
    gdat->Draw("ap");
    gpred->Sort();
    gpred->SetLineColor(kRed);
    gpred->Draw("l same");
    gdat_test->SetMarkerStyle(kFullDotMedium);
    gdat_test->SetMarkerColor(kBlue);
    gdat_test->Draw("p same");
  }

  c->cd(0);
}

Eigen::MatrixXd additive_model(const Eigen::MatrixXd& xs,
                               const Eigen::VectorXd& ys,
                               const Eigen::MatrixXd& xs_test,
                               const Eigen::VectorXd& ys_test)
{
  const unsigned int Npt = ys.size();
  const unsigned int Ntest = ys_test.size();
  assert(xs.rows() == int(Npt));
  const unsigned int Nvar = xs.cols();

  std::cout << "Solving model for " << Npt << " universes described by " << Nvar << " vars" << std::endl;

  Eigen::MatrixXd preds = Eigen::MatrixXd::Zero(Npt, Ndim);
  Eigen::MatrixXd preds_test = Eigen::MatrixXd::Zero(Ntest, Ndim);

  double old_mse = calc_mse(preds, ys);
  std::cout << "MSE " << old_mse << " (" << sqrt(old_mse) << ")" << std::endl;

  plot_preds(xs, ys, preds, xs_test, ys_test, preds_test);

  gPad->Print("preds_anim.pdf[");

  gPad->Print("preds_orig.pdf");
  gPad->Print("preds_anim.pdf");
  /*
  for(int pass = 0; pass < 100; ++pass){
    for(unsigned int ivar = 0; ivar < Nvar; ++ivar){
      // Residual
      const Eigen::VectorXd dy = ys - (total_prediction(preds) - preds[ivar]);
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

  // Start with the same basis as the regular variables
  Eigen::MatrixXd betas = Eigen::MatrixXd::Identity(Nvar, Ndim);

  /* // random initialization
  for(unsigned int i = 0; i < Nvar; ++i){
    for(unsigned int j = 0; j < Nvar; ++j){
      betas[i][j] = gRandom->Gaus();
    }
    betas[i] *= 1/sqrt(betas[i].squaredNorm()); // unit vector
  }
  */

  int pass = 0;
  //  while(true){
  for(int k = 0; k < 10000; ++k){
    for(unsigned int idim = 0; idim < Ndim; ++idim){
      // Residual
      const Eigen::VectorXd dy = ys - (total_prediction(preds) - preds.col(idim));
      Eigen::VectorXd beta = betas.col(idim);
      Eigen::VectorXd preds_test_col = preds_test.col(idim);
      preds.col(idim) = local_linear_update_basis(beta, xs, dy, xs_test, preds_test_col);
      betas.col(idim) = beta; // TODO clunky, wanted to update in place
      preds_test.col(idim) = preds_test_col;
    } // end for ivar

    if(k <= 100 || (k <= 1000 && k%10 == 9) || k%100 == 99){
      plot_preds(projectMulti(betas, xs), ys, preds,
                 projectMulti(betas, xs_test), ys_test, preds_test);
      gPad->Print(TString::Format("preds_%d.pdf", pass).Data());
      gPad->Print("preds_anim.pdf");
    }

    const double mse = calc_mse(preds, ys);
    std::cout << pass << ": MSE " << mse << " (" << sqrt(mse) << ")"
              << " test: " << sqrt(calc_mse(preds_test, ys_test))
              << std::endl;

    //    if(mse >= old_mse) break; // convergence
    if(mse >= old_mse){
      std::cout << "WOULD HAVE BROKEN HERE" << std::endl;
    }
    old_mse = mse;

    ++pass;
  } // end for pass

  gPad->Print("preds_anim.pdf]");

  //  for(const Eigen::VectorXd& beta: betas) beta.Print();

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

  Eigen::MatrixXd xs(90, systs.size());
  Eigen::MatrixXd xs_test(10, systs.size());

  for(unsigned int i = 0; i < systs.size(); ++i){
    TVectorD* v = (TVectorD*)fin->Get(TString::Format("xs_%d", i).Data());
    for(int j = 0; j < 100; ++j){
      if(j < 90)
        xs(j, i) = (*v)[j];
      if(j >= 90)
        xs_test(j-90, i) = (*v)[j];
    }
  }

  Eigen::VectorXd ys(90);
  Eigen::VectorXd ys_test(10);

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

    if(i < 90)
      ys[i] = log(y);
    if(i >= 90)
      ys_test[i-90] = log(y);
  }

  /*
  TH1* h = snom.ToTH1(pot);
  h->Draw("hist");
  h->GetYaxis()->SetRangeUser(0, 1.3*h->GetMaximum());
  for(const Spectrum& s: multiverse) s.ToTH1(pot, kGray)->Draw("hist same");
  h->Draw("hist same");
  gPad->Print("multiverse.pdf");
  */

  additive_model(xs, ys, xs_test, ys_test);
}
