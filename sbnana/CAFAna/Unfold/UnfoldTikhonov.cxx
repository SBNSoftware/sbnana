#include "sbnana/CAFAna/Unfold/UnfoldTikhonov.h"

//#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/HistCache.h"

#include "TH2.h"

#include <Eigen/Dense>

#include <cassert>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  Spectrum UnfoldTikhonov(const Spectrum& reco,
                          const ReweightableSpectrum& recoVsTrue,
                          double regStrength)
  {
    assert(regStrength >= 0);

    if(recoVsTrue.NDimensions() != 1){
      std::cout << "UnfoldTikhonov: must use 1D reco axis. This is not a fundamental limitation. The code could be improved relatively easily." << std::endl;
      abort();
    }

    DontAddDirectory guard;

    // This is where that requirement comes in
    const unsigned int Nreco = recoVsTrue.GetBinnings()[0].NBins();
    const unsigned int Ntrue = recoVsTrue.GetBinnings()[1].NBins();

    // Smearing matrix from true to reco
    Eigen::MatrixXd M(Nreco, Ntrue);

    TH2D* m = recoVsTrue.ToTH2(1);
    for(unsigned int j = 0; j < Ntrue; ++j){
      double tot = 0;
      for(unsigned int i = 0; i < Nreco; ++i){
        const double mij = m->GetBinContent(i+1, j+1);
        M(i, j) = mij; 
        tot += mij;
      }
      // Normalize each column of true energy. One true event smears to one
      // reco event total.
      if(tot != 0) for(unsigned int i = 0; i < Nreco; ++i) M(i, j) /= tot;
    } // end for j
    delete m;

    // Data vector
    TH1D* hy = reco.ToTH1(reco.POT());
    Eigen::VectorXd y(Nreco);
    for(unsigned int i = 0; i < Nreco; ++i) y[i] = hy->GetBinContent(i+1);
    HistCache::Delete(hy);

    // Matrix penalizing true distributions with large second derivative. Would
    // also need an update here for multi-dimensional distribution. I think you
    // would add an extra P term to the equations for each dimension.
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(Ntrue, Ntrue);
    for(unsigned int i = 1; i+1 < Ntrue; ++i){
      P(i, i) = -2;
      P(i, i-1) = P(i, i+1) = 1;
    }

    // We are trying to minimize a chisq
    //
    // chisq = (M*x-y)^2 + lambda*(P*x)^2
    //
    // to solve for the best true distribution x.
    //
    // This results in needing to solve
    //
    // (M^T*M + lambda*P^T*P)*x = M^T*y

    const Eigen::MatrixXd lhs = M.transpose()*M + regStrength*P.transpose()*P;
    const Eigen::VectorXd rhs = M.transpose()*y;

    // https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(lhs);
    const Eigen::VectorXd ret = dec.solve(rhs);


    // To get the correct binning
    std::unique_ptr<TH1D> hret(recoVsTrue.WeightingVariable().ToTH1(reco.POT()));
    hret->Reset();
    for(unsigned int i = 0; i < Ntrue; ++i) hret->SetBinContent(i+1, ret[i]);

    // TODO in principle these should be the true labels and bins. Will be
    // easier with cafanacore
    return Spectrum(std::move(hret), {""}, {recoVsTrue.GetBinnings()[1]},
                    reco.POT(), reco.Livetime());
  }
}
