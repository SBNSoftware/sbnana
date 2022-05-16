#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

#include "CRUMBS_helper.h"
#include "TCanvas.h"
#include "TH1D.h"

void CRUMBS_demo()
{
  const std::string inputFile = "./example.flat.caf.root";

  SpectrumLoader loader(inputFile);

  const double gPOT = 6.6e20;
  const Binning binsCRUMBSScore = Binning::Simple(20,-1,1);
  const Binning binsSliceNTracks = Binning::Simple(10,0,10);

  Spectrum specCRUMBSScore = Spectrum("CRUMBS Score", binsCRUMBSScore, loader, kCRUMBSScore, kNoSpillCut, kNonUnambiguousSlice);
  Spectrum specBestCRUMBSScore = Spectrum("Best CRUMBS Score", binsCRUMBSScore, loader, kBestCRUMBSScore, kHasNonUnambiguousSlice);
  Spectrum specBestCRUMBSSliceNTracks = Spectrum("# Tracks", binsSliceNTracks, loader, kBestCRUMBSSliceNTracks, kHasNonUnambiguousSlice && kBestCRUMBSSliceCut);

  loader.Go();

  TCanvas *canCRUMBSScore = new TCanvas("canCRUMBSScore","canCRUMBSScore");
  canCRUMBSScore->cd();

  TH1D* histCRUMBSScore = specCRUMBSScore.ToTH1(gPOT);
  histCRUMBSScore->SetLineWidth(4);
  histCRUMBSScore->SetLineColor(kMagenta+2);
  histCRUMBSScore->Draw("hist");

  TCanvas *canBestCRUMBSScore = new TCanvas("canBestCRUMBSScore","canBestCRUMBSScore");
  canBestCRUMBSScore->cd();

  TH1D* histBestCRUMBSScore = specBestCRUMBSScore.ToTH1(gPOT);
  histBestCRUMBSScore->SetLineWidth(4);
  histBestCRUMBSScore->SetLineColor(kMagenta+2);
  histBestCRUMBSScore->Draw("hist");

  TCanvas *canBestCRUMBSSliceNTracks = new TCanvas("canBestCRUMBSSliceNTracks","canBestCRUMBSSliceNTracks");
  canBestCRUMBSSliceNTracks->cd();

  TH1D* histBestCRUMBSSliceNTracks = specBestCRUMBSSliceNTracks.ToTH1(gPOT);
  histBestCRUMBSSliceNTracks->SetLineWidth(4);
  histBestCRUMBSSliceNTracks->SetLineColor(kMagenta+2);
  histBestCRUMBSSliceNTracks->Draw("hist");
}
  

  
