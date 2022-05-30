#include "cafanacore/Spectrum.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "sbnana/CAFAna/Core/Cut.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  bool IsCCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg);

  bool IsNCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg);

  class FluxTimesNuclei: public Spectrum
  {
  public:
    /// pdg PDG code for neutrino, -14,-12,+12,14
    FluxTimesNuclei(INuTruthSource& src, const Binning& bins,
                    const NuTruthCut& fidvol, int pdg, const NuTruthWeight& wgt = kNuTruthUnweighted);

    TH1D* ToTH1(double pot,
                Color_t col = kBlack,
                Style_t style = kSolid,
                EBinType bintype = kBinContent);

    void SaveTo(TDirectory* dir, const std::string& name) const;

    static std::unique_ptr<FluxTimesNuclei> LoadFrom(TDirectory* dir, const std::string& name);

    /// Convert an \ref FluxTimesNuclei into a \ref Spectrum where every bin is the integral of the
    /// \ref FluxTimesNuclei, useful for dividing out flux in cross section measurement
    Spectrum MakeTotalFlux(const LabelsAndBins& ax) const;

  protected:
    /// Helper for LoadFrom()
    FluxTimesNuclei(const std::unique_ptr<Spectrum> spec, const int pdg);

    int fPdg;
  };

  class EnsembleFluxTimesNuclei: public EnsembleSpectrum
  {
  public:
    /// pdg PDG code for neutrino, -14,-12,+12,14
    EnsembleFluxTimesNuclei(INuTruthEnsembleSource& src, const Binning& bins,
                            const NuTruthCut& fidvol, int pdg, const NuTruthWeight& wgt = kNuTruthUnweighted);

    TH1D* ToTH1(double pot,
                Color_t col = kBlack,
                Style_t style = kSolid,
                EBinType bintype = kBinContent);

    void SaveTo(TDirectory* dir, const std::string& name) const;

    static std::unique_ptr<EnsembleFluxTimesNuclei> LoadFrom(TDirectory* dir, const std::string& name);

    /// Convert an \ref EnsembleFluxTimesNuclei into a \ref EnsembleSpectrum where every bin within
    /// a given universe is the integral of the \ref EnsembleFluxTimesNuclei, useful for dividing
    /// out flux in cross section measurement
    EnsembleSpectrum MakeTotalFlux(const LabelsAndBins& ax) const;

  protected:
    /// Helper for LoadFrom()
    EnsembleFluxTimesNuclei(const std::unique_ptr<EnsembleSpectrum> spec, const int pdg);

    int fPdg;
  };


}
