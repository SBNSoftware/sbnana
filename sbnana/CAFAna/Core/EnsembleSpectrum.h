#pragma once

#include "CAFAna/Core/Spectrum.h"

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/IRecordSource.h"
#include "sbnana/CAFAna/Core/Weight.h"

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
    EnsembleSpectrum(ISliceSource& src,
                     const HistAxis& axis,
                     const std::vector<SystShifts>& univ_shifts);

    EnsembleSpectrum(ISliceSource& src,
                     const HistAxis& axis,
                     const std::vector<Weight>& univ_weis);

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
    TGraphAsymmErrors* ErrorBand(double exposure,
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

    void SaveTo(TDirectory* dir, const std::string& name) const;
    static std::unique_ptr<EnsembleSpectrum> LoadFrom(TDirectory* dir,
                                                      const std::string& name);

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
