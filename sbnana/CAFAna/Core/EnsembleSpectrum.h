#pragma once

#include "sbnana/CAFAna/Core/Spectrum.h"

#include <cassert>

class TGraphAsymmErrors;

namespace ana
{
  class EnsembleRatio;

  // TODO Multiverse class encapsulating a vector of shifts/weights and an ID
  // number?

  class EnsembleSpectrum
  {
  public:
    EnsembleSpectrum(SpectrumLoaderBase& loader,
                     const HistAxis& axis,
                     const SpillCut& spillcut,
                     const Cut& cut,
                     const std::vector<SystShifts>& univ_shifts,
                     const Var& cv_wei = kUnweighted);

    EnsembleSpectrum(SpectrumLoaderBase& loader,
                     const HistAxis& axis,
                     const SpillCut& spillcut,
                     const Cut& cut,
                     const std::vector<Var>& univ_weis,
                     const Var& cv_wei = kUnweighted);

    EnsembleSpectrum(SpectrumLoaderBase& loader,
                     const HistAxis& xAxis,
                     const HistAxis& yAxis,
                     const SpillCut& spillcut,
                     const Cut& cut,
                     const std::vector<SystShifts>& univ_shifts,
                     const Var& cv_wei = kUnweighted);

    EnsembleSpectrum(SpectrumLoaderBase& loader,
                     const HistAxis& xAxis,
                     const HistAxis& yAxis,
                     const SpillCut& spillcut,
                     const Cut& cut,
                     const std::vector<Var>& univ_weis,
                     const Var& cv_wei = kUnweighted);

    Spectrum Nominal() const {return fNom;}
    unsigned int NUniverses() const {return fUnivs.size();}
    Spectrum Universe(unsigned int i) const
    {
      assert(i < fUnivs.size());
      return fUnivs[i];
    }

    double POT() const {return fNom.POT();}

    double Livetime() const {return fNom.Livetime();}

    /// Result can be painted prettily with \ref DrawErrorBand
    TGraphAsymmErrors* ErrorBand(double z,
                                 double exposure,
                                 EExposureType expotype = kPOT,
                                 EBinType bintype = kBinContent) const;

    void Scale(double c);

    EnsembleSpectrum& operator+=(const EnsembleSpectrum& rhs);
    EnsembleSpectrum operator+(const EnsembleSpectrum& rhs) const;

    EnsembleSpectrum& operator-=(const EnsembleSpectrum& rhs);
    EnsembleSpectrum operator-(const EnsembleSpectrum& rhs) const;

    EnsembleSpectrum& operator*=(const EnsembleRatio& rhs);
    EnsembleSpectrum operator*(const EnsembleRatio& rhs) const;

    EnsembleSpectrum& operator/=(const EnsembleRatio& rhs);
    EnsembleSpectrum operator/(const EnsembleRatio& rhs) const;

    void SaveTo(TDirectory* dir) const;
    static std::unique_ptr<EnsembleSpectrum> LoadFrom(TDirectory* dir);

    unsigned int NDimensions() const{return fNom.NDimensions();}
    std::vector<std::string> GetLabels() const {return fNom.GetLabels();}
    std::vector<Binning> GetBinnings() const {return fNom.GetBinnings();}

  protected:
    EnsembleSpectrum(const Spectrum& nom) : fNom(nom) {}

    Spectrum fNom;
    std::vector<Spectrum> fUnivs;
  };

  // Commutative
  inline EnsembleSpectrum operator*(const EnsembleRatio& lhs, const EnsembleSpectrum& rhs){return rhs*lhs;}
  inline EnsembleSpectrum operator/(const EnsembleRatio& lhs, const EnsembleSpectrum& rhs){return rhs/lhs;}

  void DrawErrorBand(TH1* nom,
                     TGraphAsymmErrors* band,
                     int bandCol = -1,
                     double alpha = 1);
}
