#include "sbnana/CAFAna/XSec/Flux.h"

#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Weight.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1.h"
#include "TVectorD.h"

namespace ana
{
  //----------------------------------------------------------------------
  bool IsCCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg)
  {
    // It doesn't matter exactly what process we choose, so long as it
    // corresponds to one of the modes GENIE considers as a top-level
    // cross-section. Here we go for CC QE on Argon.
    return nu->pdg == pdg && nu->iscc &&
      nu->genie_mode == caf::kQE && nu->genie_inttype == caf::kCCQE &&
      nu->targetPDG == 1000180400 /* Argon 40 */ &&
      !nu->ischarm;
  }

  //----------------------------------------------------------------------
  bool IsNCQEOnArgon(const caf::SRTrueInteractionProxy* nu, int pdg)
  {
    // The issue with the CC interaction is there is a threshold around the
    // muon mass. This process works for lower neutrino energies too
    return nu->pdg == pdg && !nu->iscc &&
      nu->genie_mode == caf::kQE && nu->genie_inttype == caf::kNCQE &&
      nu->targetPDG == 1000180400 /* Argon 40 */ &&
      !nu->ischarm &&
      nu->hitnuc == 2112;
  }

  //----------------------------------------------------------------------
  const NuTruthWeight kInvXSec([](const caf::SRTrueInteractionProxy* nu)
                             {
                               // GENIE uses GeV internally. We ultimately want
                               // a flux in m^-2
                               const double GeV2perm2 = 2.56819e31;

                               return GeV2perm2/nu->xsec;
                             });

  //----------------------------------------------------------------------
  NuTruthCut IsNCQEOnArgonCut(int pdg)
  {
    return NuTruthCut([pdg](const caf::SRTrueInteractionProxy* nu)
                      {
                        return IsNCQEOnArgon(nu, pdg);
                      });
  }

  //----------------------------------------------------------------------
  FluxTimesNuclei::FluxTimesNuclei(INuTruthSource& src,
                                   const Binning& bins,
                                   const NuTruthCut& fidvol,
                                   int pdg,
                                   const NuTruthWeight& wgt)
    : Spectrum(src[IsNCQEOnArgonCut(pdg) && fidvol].Weighted(wgt *kInvXSec),
               NuTruthHistAxis("True neutrino energy (GeV)",
                               bins,
                               SIMPLENUTRUTHVAR(E))),
      fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  TH1D* FluxTimesNuclei::ToTH1(double pot,
                               Color_t col,
                               Style_t style,
                               EBinType bintype)
  {
    TH1D* ret = Spectrum::ToTH1(pot, col, style, kPOT, bintype);

    std::string ytitle = "Flux #times nuclei (";
    switch(fPdg){
    case +12: ytitle += "#nu_e"; break;
    case -12: ytitle += "#bar{#nu}_e"; break;
    case +14: ytitle += "#nu_{#mu}"; break;
    case -14: ytitle += "#bar{#nu}_{#mu}"; break;
    }

    ytitle += "/ m^{2}";
    if(bintype == kBinDensity) ytitle += " / GeV";
    ytitle += ")";
    ret->GetYaxis()->SetTitle(ytitle.c_str());
    return ret;
  }

  //----------------------------------------------------------------------
  void FluxTimesNuclei::SaveTo(TDirectory* dir, const std::string& name) const
  {
    TDirectory* tmp = gDirectory;

    dir = dir->mkdir(name.c_str()); // switch to subdir
    dir->cd();

    TObjString("FluxTimesNuclei").Write("type");

    TVectorD pdg(1);
    pdg[0] = fPdg;
    pdg.Write("pdg");

    Spectrum::SaveTo(dir, "spectrum");

    dir->Write();
    delete dir;

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<FluxTimesNuclei> FluxTimesNuclei::LoadFrom(TDirectory* dir, const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str()); // switch to subdir
    assert(dir);

    DontAddDirectory guard;

    std::unique_ptr<TObjString> tag((TObjString*)dir->Get("type"));
    assert(tag);
    assert(tag->GetString() == "FluxTimesNuclei");

    std::unique_ptr<TVectorD> pdg((TVectorD*)dir->Get("pdg"));
    assert(pdg && pdg->GetNrows()==1);

    return std::unique_ptr<FluxTimesNuclei>(new FluxTimesNuclei(Spectrum::LoadFrom(dir, "spectrum"), (*pdg)[0]));
  }

  //----------------------------------------------------------------------
  FluxTimesNuclei::FluxTimesNuclei(const std::unique_ptr<Spectrum> spec, const int pdg)
    : Spectrum(*spec), fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  FluxTimesNuclei::FluxTimesNuclei(const Eigen::ArrayXd&& hist,
                                   const LabelsAndBins& axis,
                                   double pot,
                                   double livetime,
                                   int pdg)
    : Spectrum(hist, axis, pot, livetime), fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  FluxTimesNuclei FluxTimesNuclei::MakeTotalFlux(const LabelsAndBins& axis) const
  {
    const unsigned int nbins = axis.GetBins1D().NBins()+2;

    const Eigen::ArrayXd hist = Eigen::ArrayXd::Constant(nbins, Integral(fPOT));

    return FluxTimesNuclei(std::move(hist), axis, fPOT, fLivetime, fPdg);
  }

  //----------------------------------------------------------------------
  EnsembleFluxTimesNuclei::EnsembleFluxTimesNuclei(INuTruthEnsembleSource& src,
                                                   const Binning& bins,
                                                   const NuTruthCut& fidvol,
                                                   int pdg,
                                                   const NuTruthWeight& wgt)
    : EnsembleSpectrum(src[IsNCQEOnArgonCut(pdg) && fidvol].Weighted(wgt * kInvXSec),
                       NuTruthHistAxis("True neutrino energy (GeV)",
                                       bins,
                                       SIMPLENUTRUTHVAR(E))),
      fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  TH1D* EnsembleFluxTimesNuclei::ToTH1(double pot,
                                       Color_t col,
                                       Style_t style,
                                       EBinType bintype)
  {
    TH1D* ret = EnsembleSpectrum::Nominal().ToTH1(pot, col, style, kPOT, bintype);

    std::string ytitle = "Flux #times nuclei (";
    switch(fPdg){
    case +12: ytitle += "#nu_e"; break;
    case -12: ytitle += "#bar{#nu}_e"; break;
    case +14: ytitle += "#nu_{#mu}"; break;
    case -14: ytitle += "#bar{#nu}_{#mu}"; break;
    }

    ytitle += "/ m^{2}";
    if(bintype == kBinDensity) ytitle += " / GeV";
    ytitle += ")";
    ret->GetYaxis()->SetTitle(ytitle.c_str());
    return ret;
  }

  //----------------------------------------------------------------------
  void EnsembleFluxTimesNuclei::SaveTo(TDirectory* dir, const std::string& name) const
  {
    TDirectory* tmp = gDirectory;

    dir = dir->mkdir(name.c_str()); // switch to subdir
    dir->cd();

    TObjString("EnsembleFluxTimesNuclei").Write("type");

    TVectorD pdg(1);
    pdg[0] = fPdg;
    pdg.Write("pdg");

    EnsembleSpectrum::SaveTo(dir, "ensemblespectrum");

    dir->Write();
    delete dir;

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<EnsembleFluxTimesNuclei> EnsembleFluxTimesNuclei::LoadFrom(TDirectory* dir, const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str()); // switch to subdir
    assert(dir);

    DontAddDirectory guard;

    std::unique_ptr<TObjString> tag((TObjString*)dir->Get("type"));
    assert(tag);
    assert(tag->GetString() == "EnsembleFluxTimesNuclei");

    std::unique_ptr<TVectorD> pdg((TVectorD*)dir->Get("pdg"));
    assert(pdg && pdg->GetNrows()==1);

    return std::unique_ptr<EnsembleFluxTimesNuclei>(new EnsembleFluxTimesNuclei(EnsembleSpectrum::LoadFrom(dir, "ensemblespectrum"), (*pdg)[0]));
  }

  //----------------------------------------------------------------------
  EnsembleFluxTimesNuclei::EnsembleFluxTimesNuclei(const std::unique_ptr<EnsembleSpectrum> spec, const int pdg)
    : EnsembleSpectrum(*spec), fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  EnsembleFluxTimesNuclei::EnsembleFluxTimesNuclei(const FitMultiverse* multiverse,
                                                   const Hist&& hist,
                                                   double pot,
                                                   double livetime,
                                                   const LabelsAndBins& axis,
                                                   int pdg)
    : EnsembleSpectrum(multiverse, std::move(hist), pot, livetime, axis), fPdg(pdg)
  {
  }

  //----------------------------------------------------------------------
  EnsembleFluxTimesNuclei EnsembleFluxTimesNuclei::MakeTotalFlux(const LabelsAndBins& axis) const
  {
    const unsigned int nuniv = NUniverses();
    const unsigned int nbins = axis.GetBins1D().NBins()+2;

    Hist hist = Hist::Zero(nbins * nuniv);

    for(unsigned int univIdx = 0; univIdx < nuniv; ++univIdx){
      const double univIntegral(Universe(univIdx).Integral(fPOT));
      for(unsigned int bin = 0; bin < nbins; ++bin)
        hist.Fill(nbins * univIdx + bin, univIntegral);
    }

    return EnsembleFluxTimesNuclei(fMultiverse, std::move(hist), fPOT, fLivetime, axis, fPdg);
  }
}
