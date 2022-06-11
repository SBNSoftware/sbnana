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
    FluxTimesNuclei MakeTotalFlux(const LabelsAndBins& axis) const;

    int PDG() const {return fPdg;};

  protected:
    /// Helper for LoadFrom()
    FluxTimesNuclei(const Spectrum&& spec, const int pdg);
    FluxTimesNuclei(const Eigen::ArrayXd&& hist,
                    const LabelsAndBins& axis,
                    double pot,
                    double livetime,
                    int pdg);

    int fPdg;

    friend class EnsembleFluxTimesNuclei;
  };

  class EnsembleFluxTimesNuclei: public EnsembleSpectrum
  {
  public:
    /// pdg PDG code for neutrino, -14,-12,+12,14
    EnsembleFluxTimesNuclei(INuTruthEnsembleSource& src, const Binning& bins,
                            const NuTruthCut& fidvol, int pdg, const NuTruthWeight& wgt = kNuTruthUnweighted);

    /// \brief Creates an ensemble flux times nuceli from a nominal input \ref FluxTimesNuclei
    //         which is replicated nUniverse times from the multiverse which it adopts.
    //         Note that this is a temporary workaround for now
    static EnsembleFluxTimesNuclei ReplicatedNominal(const FluxTimesNuclei& spec, const FitMultiverse* multiverse);

    TH1D* ToTH1(double pot,
                Color_t col = kBlack,
                Style_t style = kSolid,
                EBinType bintype = kBinContent);

    void SaveTo(TDirectory* dir, const std::string& name) const;

    static std::unique_ptr<EnsembleFluxTimesNuclei> LoadFrom(TDirectory* dir, const std::string& name);
    FluxTimesNuclei Nominal() const {return Universe(0);}
    FluxTimesNuclei Universe(unsigned int i) const{return FluxTimesNuclei(EnsembleSpectrum::Universe(i), fPdg);};

    /// Convert an \ref EnsembleFluxTimesNuclei into a \ref EnsembleSpectrum where every bin within
    /// a given universe is the integral of the \ref EnsembleFluxTimesNuclei, useful for dividing
    /// out flux in cross section measurement
    EnsembleFluxTimesNuclei MakeTotalFlux(const LabelsAndBins& axis) const;

    int PDG() const {return fPdg;};

  protected:
    /// Helper for LoadFrom()
    EnsembleFluxTimesNuclei(const EnsembleSpectrum&& spec, const int pdg);
    EnsembleFluxTimesNuclei(const FitMultiverse* multiverse,
                            const Hist&& hist,
                            double pot,
                            double livetime,
                            const LabelsAndBins& axis,
                            int pdg);

    int fPdg;
  };
}
