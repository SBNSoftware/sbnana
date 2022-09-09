#include "sbnana/CAFAna/Unfold/UnfoldSVD.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "TSVDUnfold.h"

#include "TH2.h"

namespace ana
{
  //----------------------------------------------------------------------
  Spectrum UnfoldSVD(const Spectrum& reco,
                     const ReweightableSpectrum& recoVsTrue,
                     unsigned int reg)
  {
    const double pot = reco.POT();

    DontAddDirectory guard;

    std::unique_ptr<TH1D> hReco(reco.ToTH1(pot));
    std::unique_ptr<TH1D> hMCReco(recoVsTrue.UnWeighted().ToTH1(pot));
    std::unique_ptr<TH1D> hMCTrue(recoVsTrue.WeightingVariable().ToTH1(pot));
    std::unique_ptr<TH2D> hRT(recoVsTrue.ToTH2(pot));

    TSVDUnfold uf(hReco.get(), hMCReco.get(), hMCTrue.get(), hRT.get());

    std::unique_ptr<TH1D> h_unf(uf.Unfold(reg));

    // Enforce matching normalization
    h_unf->Scale(reco.Integral(pot)/h_unf->Integral(0, -1));

    // TODO in principle these should be the true labels and bins. Will be
    // easier with cafanacore
    Eigen::ArrayXd arr = Eigen::Map<Eigen::ArrayXd>(h_unf->GetArray(),
                                                    h_unf->GetNbinsX()+2);
    return Spectrum(arr,
                    LabelsAndBins(reco.GetLabels(), recoVsTrue.GetTrueBinnings()),
                    pot, reco.Livetime());
  }
}
