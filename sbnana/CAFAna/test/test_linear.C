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


void linear_model(const Eigen::MatrixXd& xs,
                  const Eigen::VectorXd& ys,
                  const Eigen::MatrixXd& xs_test,
                  const Eigen::VectorXd& ys_test)
{
  std::cout << "Orig MSE: " << ys.squaredNorm() / ys.size() << " ("
            << sqrt(ys.squaredNorm() / ys.size()) << ")" << std::endl;
  std::cout << "test: " << ys_test.squaredNorm() / ys_test.size() << "("
            << sqrt(ys_test.squaredNorm() / ys_test.size()) << ")" << std::endl;

  // A is vec of coefficients
  // xs * A = ys
  // minimize (xs*A - ys)^2
  // d/dA -> 2*(xs*A - ys)*xs = 0
  // xs*xsT*A = ys*xs

  const Eigen::VectorXd A = (xs.transpose()*xs).ldlt().solve(xs.transpose()*ys);  

  std::cout << A << std::endl;

  const Eigen::VectorXd preds = xs * A;

  const double mse = (preds-ys).squaredNorm() / ys.size();

  std::cout << "MSE = " << mse << " (" << sqrt(mse) << ")" << std::endl;

  const Eigen::VectorXd preds_test = xs_test * A;
  const double mse_test = (preds_test-ys_test).squaredNorm() / ys_test.size();

  std::cout << "test MSE = " << mse_test << " (" << sqrt(mse_test) << ")" << std::endl;
}


void test_linear()
{
  TFile* fin = new TFile("additive_state.root");
  Spectrum snom = *LoadFrom<Spectrum>(fin->GetDirectory("nominal"));
  std::vector<Spectrum> multiverse;
  for(int i = 0; i < 100; ++i) multiverse.push_back(*LoadFrom<Spectrum>(fin->GetDirectory(TString::Format("multiverse_%d", i).Data())));

  Eigen::MatrixXd xs = Eigen::MatrixXd::Zero(90, systs.size());
  Eigen::MatrixXd xs_test = Eigen::MatrixXd::Zero(10, systs.size());

  for(unsigned int i = 0; i < systs.size(); ++i){
    //    if(i == 1 || i == 2) continue; // smallest coeffs
    //    if(i == 0 || i == 3 || i == 13 || i == 14 || i == 16 || i == 17 || i == 19 || i == 22 || i == 25) continue; // next smallest
    //    if(i == 5 || i == 8 || i == 12 || i == 15 || i == 18 || i == 23 || i == 24) continue;
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

    if(i < 90)
      ys[i] = log(y);
    if(i >= 90)
      ys_test[i-90] = log(y);
  }

  linear_model(xs, ys, xs_test, ys_test);
}
