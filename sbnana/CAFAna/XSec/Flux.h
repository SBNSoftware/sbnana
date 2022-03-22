#include "sbnana/CAFAna/Core/Spectrum.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  bool IsCCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg);

  bool IsNCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg);

  class FluxTimesNuclei: public Spectrum
  {
  public:
    /// pdg PDG code for neutrino, -14,-12,+12,14
    FluxTimesNuclei(SpectrumLoaderBase& loader, const Binning& bins,
                    const Cut& fidvol, int pdg);

    TH1D* ToTH1(double exposure,
                Color_t col = kBlack,
                Style_t style = kSolid,
                EBinType bintype = kBinContent);
  protected:
    int fPdg;
  };
}
