#pragma once

#include "sbnana/CAFAna/Analysis/Surface.h"

#include <vector>

namespace ana
{
  class MedianSurface: public Surface
  {
  public:
    MedianSurface(const std::vector<Surface>& throws);

    void DrawEnsemble(TH2* fc, Color_t color = kGray);
    void DrawBand(TH2* fc);

    void SaveTo(TDirectory* dir, const std::string& name) const;
    static std::unique_ptr<MedianSurface> LoadFrom(TDirectory* dir,
                                                   const std::string& name);
  protected:
    std::vector<Surface> fThrows;

    TH2F *fHistUp1, *fHistDn1, *fHistUp2, *fHistDn2;
  };
}
