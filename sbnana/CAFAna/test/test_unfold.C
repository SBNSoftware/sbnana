#include "sbnana/CAFAna/Unfold/UnfoldIterative.h"
#include "sbnana/CAFAna/Unfold/UnfoldSVD.h"
#include "sbnana/CAFAna/Unfold/UnfoldTikhonov.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TRandom3.h"

using namespace ana;

const Binning recobins = Binning::Simple(40, 0, 10);
// Must be the same as recobins if you want SVD unfolding to work
const Binning truebins = Binning::Simple(80, 0, 20);

ReweightableSpectrum GetRecoVsTrue()
{
  // Mock up a simple ReweightableSpectrum with a gaussian true distribution
  // and a reco distribution that is shifted and scaled. In a real analysis you
  // would fill your ReweightableSpectrum using SpectrumLoader in the usual
  // way.

  TH2* hrvt = new TH2F("", "Smearing matrix;Reco;True",
                       recobins.NBins(), recobins.Min(), recobins.Max(),
                       truebins.NBins(), truebins.Min(), truebins.Max());

  for(int i = 0; i < 1e6; ++i){
    const double t = gRandom->Gaus(5, 2);
    const double r = t*gRandom->Gaus(.7, .2);
    hrvt->Fill(r, t);
  }

  return ReweightableSpectrum(kUnweighted, hrvt,
                              {"Reco"},
                              {recobins, truebins},
                              500, 0);

}

void DoIterative(const Spectrum& reco, const Spectrum& truth,
                 const ReweightableSpectrum& rvt)
{
  // Iterative unfolding with decreasing regularization strength

  const Spectrum it1 = UnfoldIterative(reco, rvt, 1);
// const Spectrum it3 = UnfoldIterative(reco, rvt, 3);
// const  Spectrum it5 = UnfoldIterative(reco, rvt, 5);
  const Spectrum it10 = UnfoldIterative(reco, rvt, 10);

  TH1* hreco = reco.ToTH1(1);
  hreco->SetTitle("Iterative unfolding");
  hreco->Draw("ep");
  truth.ToTH1(1, kRed)->Draw("hist same");
  it1.ToTH1(1, kBlue)->Draw("hist same");
//  it3.ToTH1(1, kBlue)->Draw("hist same");
//  it5.ToTH1(1, kBlue)->Draw("hist same");
  it10.ToTH1(1, kBlue-7)->Draw("hist same");
}

void DoSVD(const Spectrum& reco, const Spectrum& truth,
           const ReweightableSpectrum& rvt)
{
  // Singular value decomposition unfolding with decreaseing regularization

  const Spectrum svd2 = UnfoldSVD(reco, rvt, 2);
  const Spectrum svd4 = UnfoldSVD(reco, rvt, 4);
  const Spectrum svd10 = UnfoldSVD(reco, rvt, 10);

  TH1* hreco = reco.ToTH1(1);
  hreco->SetTitle("SVD unfolding");
  hreco->Draw("ep");
  truth.ToTH1(1, kRed)->Draw("hist same");
  svd2.ToTH1(1, kBlue)->Draw("hist same");
  svd4.ToTH1(1, kBlue-3)->Draw("hist same");
  svd10.ToTH1(1, kBlue-7)->Draw("hist same");
}

void DoTikhonov(const Spectrum& reco, const Spectrum& truth,
                const ReweightableSpectrum& rvt)
{
  // Tikhonov unfolding with increasing regularization strength

  const Spectrum tikp1 = UnfoldTikhonov(reco, rvt, 1e-2);
  const Spectrum tik1 = UnfoldTikhonov(reco, rvt, 1);
  const Spectrum tik10 = UnfoldTikhonov(reco, rvt, 1e2);

  TH1* hreco = reco.ToTH1(1);
  hreco->SetTitle("Tikhonov unfolding");
  hreco->Draw("ep");
  truth.ToTH1(1, kRed)->Draw("hist same");
  tikp1.ToTH1(1, kBlue)->Draw("hist same");
  tik1.ToTH1(1, kBlue-3)->Draw("hist same");
  tik10.ToTH1(1, kBlue-7)->Draw("hist same");

  // Demonstrate that the under-regularized fit still reproduces the reco
  // spectrum correctly. Warning: this *updates rvt in-place* so you can only
  // do it once, at the end of the function, and have to change rvt to be
  // pass-by-value.

  //  rvt.ReweightToTrueSpectrum(tikp1);
  //  rvt.UnWeighted().ToTH1(1, kMagenta)->Draw("hist same");
}

void test_unfold()
{
  const ReweightableSpectrum rvt = GetRecoVsTrue();

  rvt.ToTH2(1)->Draw("colz");

  // The true distribution we are trying to retrieve by unfolding
  const Spectrum truth = rvt.WeightingVariable();
  // A stastically fluctuated reco distribution
  const Spectrum reco = rvt.UnWeighted().MockData(1);

  new TCanvas;
  DoIterative(reco, truth, rvt);

  if(recobins.NBins() == truebins.NBins()){
    new TCanvas;
    DoSVD(reco, truth, rvt);
  }
  else{
    std::cout << "Skipping SVD since reco and true bins are not matched" << std::endl;
  }

  new TCanvas;
  DoTikhonov(reco, truth, rvt);
}
