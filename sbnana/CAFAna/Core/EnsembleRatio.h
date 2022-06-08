#pragma once

#include "cafanacore/Ratio.h"

#include "TGraphAsymmErrors.h"
#include "TMatrixD.h"

#include <cassert>

namespace ana
{
  class EnsembleSpectrum;
  class FitMultiverse;

  class EnsembleRatio
  {
  public:
    friend class EnsembleSpectrum;

    EnsembleRatio(const EnsembleSpectrum& num,
                  const EnsembleSpectrum& denom,
                  bool purOrEffErrs = false);

    Ratio Nominal() const {return Universe(0);}
    unsigned int NUniverses() const;
    Ratio Universe(unsigned int i) const;

    const FitMultiverse& GetMultiverse() const {return *fMultiverse;}

    TGraphAsymmErrors* ErrorBand() const;

    /** \brief Compute bin-to-bin covariance matrix from a collection of sets of bin contents.

      \param firstBin  The first bin that should be considered (inclusive)
      \param lastBin   The last bin that should be considered (inclusive).  -1 means "last in set"

      \returns  unique_ptr to TMatrixD containing computed covariance matrix unless binSets.size() < 2,
      in which case the unique_ptr's target is nullptr.

      Note TH1D is a child class of TArrayD -- so you can pass a vector
      of TH1D* to this method.
     **/
    std::unique_ptr<TMatrixD> CalcCovMx(const int firstBin=0, const int lastBin=-1);

    /** \brief Compute bias from a collection of sets of bin contents.

      \param firstBin  The first bin that should be considered (inclusive)
      \param lastBin   The last bin that should be considered (inclusive).  -1 means "last in set"

      \returns  unique_ptr to TMatrixD containing computed bias matrix unless binSets.size() < 2,
      in which case the unique_ptr's target is nullptr.

      Note TH1D is a child class of TArrayD -- so you can pass a vector
      of TH1D* to this method.
     **/
    std::unique_ptr<TMatrixD> CalcBiasMx(const int firstBin=0, const int lastBin=-1);

    EnsembleRatio& operator*=(const EnsembleRatio& rhs);
    EnsembleRatio operator*(const EnsembleRatio& rhs) const;

    EnsembleRatio& operator/=(const EnsembleRatio& rhs);
    EnsembleRatio operator/(const EnsembleRatio& rhs) const;

  protected:
    friend class EnsembleReweightableSpectrum;
    Eigen::ArrayXd GetEigen() const {return fHist.GetEigen();}

    void CheckMultiverses(const FitMultiverse& rhs,
                          const std::string& func) const;

    const FitMultiverse* fMultiverse;
    Hist fHist;
    LabelsAndBins fAxis;
  };

  inline EnsembleRatio operator/(const EnsembleSpectrum& lhs,
                                 const EnsembleSpectrum& rhs)
  {
    return EnsembleRatio(lhs, rhs);
  }
}
