#include "sbnana/CAFAna/Unfold/UnfoldIterative.h"
#include "sbnana/CAFAna/Unfold/UnfoldSVD.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TRandom3.h"

using namespace ana;

void test_unfold()
{
  TH2* hrvt = new TH2F("", ";reco;true", 40, 0, 10, 40, 0, 10);

  for(int i = 0; i < 1e6; ++i){
    double t = gRandom->Gaus(5, 2);
    double r = t*gRandom->Gaus(.7, .2);
    hrvt->Fill(r, t);
  }

  hrvt->Draw("colz");

  Binning recobins = Binning::Simple(40, 0, 10);
  Binning truebins = recobins; // UnfoldSVD requires these to match // Binning::Simple(80, 0, 20);

  ReweightableSpectrum rvt(kUnweighted, hrvt,
                           {"Reco"},//, "True"},
                           {recobins, truebins},
                           500, 0);

  new TCanvas;

  Spectrum reco = rvt.UnWeighted().MockData(1);
  Spectrum truth = rvt.WeightingVariable();

  Spectrum it1 = UnfoldIterative(reco, rvt, 1);
//  Spectrum it3 = UnfoldIterative(reco, rvt, 3);
//  Spectrum it5 = UnfoldIterative(reco, rvt, 5);
  Spectrum it10 = UnfoldIterative(reco, rvt, 10);

  reco.ToTH1(1)->Draw("ep");
  truth.ToTH1(1, kRed)->Draw("hist same");
  it1.ToTH1(1, kBlue)->Draw("hist same");
//  it3.ToTH1(1, kBlue)->Draw("hist same");
//  it5.ToTH1(1, kBlue)->Draw("hist same");
  it10.ToTH1(1, kBlue-7)->Draw("hist same");

  new TCanvas;

  Spectrum svd2 = UnfoldSVD(reco, rvt, 2);
  Spectrum svd4 = UnfoldSVD(reco, rvt, 4);
  Spectrum svd10 = UnfoldSVD(reco, rvt, 10);
  reco.ToTH1(1)->Draw("ep");
  truth.ToTH1(1, kRed)->Draw("hist same");
  svd2.ToTH1(1, kBlue)->Draw("hist same");
  svd4.ToTH1(1, kBlue-3)->Draw("hist same");
  svd10.ToTH1(1, kBlue-7)->Draw("hist same");
}
