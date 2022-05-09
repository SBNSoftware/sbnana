#include "sbnana/CAFAna/Analysis/MedianSurface.h"

#include "sbnana/CAFAna/Core/LoadFromFile.h"

#include "TDirectory.h"
#include "TGraph.h"
#include "TH2.h"
#include "TObjString.h"
#include "TPad.h"
#include "TStyle.h"

namespace ana
{
  // --------------------------------------------------------------------------
  MedianSurface::MedianSurface(const std::vector<Surface>& throws)
    : fThrows(throws)
  {
    assert(!throws.empty());

    for(const Surface& s: throws){
      assert(s.fHist->GetNbinsX() == throws[0].fHist->GetNbinsX());
      assert(s.fHist->GetNbinsY() == throws[0].fHist->GetNbinsY());
      // TODO check min and max match at least
    }

    fLogX = throws[0].fLogX;
    fLogY = throws[0].fLogY;

    // Think about what to do with these
    fMinX = -1;
    fMinY = -1;
    // Though this is right
    fMinChi = 0;

    fHist = new TH2F(*throws[0].fHist);
    fHist->Reset();
    fHistUp1 = new TH2F(*fHist);
    fHistDn1 = new TH2F(*fHist);
    fHistUp2 = new TH2F(*fHist);
    fHistDn2 = new TH2F(*fHist);

    for(int ix = 0; ix < fHist->GetNbinsX()+2; ++ix){
      for(int iy = 0; iy < fHist->GetNbinsY()+2; ++iy){
        std::vector<float> chis;
        chis.reserve(throws.size());
        for(const Surface& s: throws){
          chis.push_back(s.fHist->GetBinContent(ix, iy));
        }
        std::sort(chis.begin(), chis.end());

        // TODO think about off-by-one errors and interpolation here
        // Median
        fHist->SetBinContent(ix, iy, chis[throws.size()/2]);
        // One sigma bounds
        const double tail1 = (1-0.6827)/2;
        fHistUp1->SetBinContent(ix, iy, chis[throws.size()*tail1]);
        fHistDn1->SetBinContent(ix, iy, chis[throws.size()*(1-tail1)]);

        // Two sigma bounds
        const double tail2 = (1-0.9545)/2;
        fHistUp2->SetBinContent(ix, iy, chis[throws.size()*tail2]);
        fHistDn2->SetBinContent(ix, iy, chis[(throws.size()-1)*(1-tail2)]);
      } // end for iy
    } // end for ix
  }

  // --------------------------------------------------------------------------
  void MedianSurface::DrawEnsemble(TH2* fc, Color_t color)
  {
    EnsureAxes();

    for(Surface& s: fThrows){
      std::vector<TGraph*> gs = s.GetGraphs(fc, -1);

      for(TGraph* g: gs){
        g->SetLineWidth(1);
        g->SetLineColor(color);
        g->Draw("l");
      }
    }

    gPad->Update();
  }

  // --------------------------------------------------------------------------
  void MedianSurface::DrawBand(TH2* fc)
  {
    EnsureAxes();

    TH2F surf1(*fHist);
    TH2F surf2(*fHist);

    for(int ix = 0; ix < fHist->GetNbinsX()+2; ++ix){
      for(int iy = 0; iy < fHist->GetNbinsY()+2; ++iy){
        const double c = fc->GetBinContent(ix, iy);

        const double u1 = fHistUp1->GetBinContent(ix, iy);
        const double d1 = fHistDn1->GetBinContent(ix, iy);
        surf1.SetBinContent(ix, iy, std::max(u1-c, c-d1));

        const double u2 = fHistUp2->GetBinContent(ix, iy);
        const double d2 = fHistDn2->GetBinContent(ix, iy);
        surf2.SetBinContent(ix, iy, std::max(u2-c, c-d2));
      }
    }

    const double level = 0;

    // I haven't been able to figure out how to draw these filled properly. cont0 does it, but can't handle countours that reach the edge of the space
    surf2.SetContour(1, &level);
    surf2.SetLineColor(kYellow);
    surf2.DrawCopy("cont3 same");

    surf1.SetContour(1, &level);
    surf1.SetLineColor(kGreen);
    surf1.SetFillColor(kGreen+2);
    surf1.DrawCopy("cont3 same");

    gPad->Update();
  }

  //----------------------------------------------------------------------
  void MedianSurface::SaveTo(TDirectory* dir, const std::string& name) const
  {
    TDirectory* tmp = gDirectory;

    dir = dir->mkdir(name.c_str()); // switch to subdir
    dir->cd();

    TObjString("MedianSurface").Write("type");

    for(unsigned int i = 0; i < fThrows.size(); ++i){
      fThrows[i].SaveTo(dir, TString::Format("surf%d", i).Data());
    }

    dir->Write();
    delete dir;

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<MedianSurface> MedianSurface::LoadFrom(TDirectory* dir,
                                                         const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str()); // switch to subdir
    assert(dir);

    DontAddDirectory guard;

    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "MedianSurface");

    std::vector<Surface> surfs;
    for(unsigned int i = 0; ; ++i){
      const std::string subname = TString::Format("surf%d", i).Data();
      TDirectory* subdir = dir->GetDirectory(subname.c_str());
      if(!subdir) break; // we got all of them
      delete subdir;
      surfs.push_back(*ana::LoadFrom<Surface>(dir, subname));
    }

    return std::make_unique<MedianSurface>(surfs);
  }

}
