#pragma once

#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include <string>

namespace ana
{
  /// \brief Make a series of spectra leaving out one cut in turn
  ///
  /// This function does *not* call Go() on the loader, allowing you to combine
  /// it with other spectra-filling if desired.
  ///
  /// \param loader  SpectrumLoader to associate spectra with
  /// \param sigcut  Definition of signal (inverse defines background)
  /// \param presel  Cuts that will always be applied
  /// \param axes    Variables that will be cut on
  /// \param cut_pos Cut position for each variable
  /// \param[out] sigs Resulting signal spectra
  /// \param[out] bkgs Resulting background spectra
  void MakeNMinusOneSpectra(SpectrumLoader& loader,
                            const Cut& sigcut,
                            const Cut& presel,
                            const std::vector<HistAxis>& axes,
                            const std::vector<double>& cut_pos,
                            std::vector<Spectrum>& sigs,
                            std::vector<Spectrum>& bkgs);

  /// \brief Search for optimum cut position given signal and background histograms
  ///
  /// \param hsig Signal distribution
  /// \param hbkg Background distribution
  /// \param[out] best_fom Figure of Merit obtained by best cut
  /// \returns Best cut position
  double FindOptimumCut(TH1* hsig, TH1* hbkg, double& best_fom);

  /// \brief Scan all cuts and update the one giving the largest FOM gain
  ///
  /// \param wildcard File name / wildcard / dataset name
  /// \param pot POT to scale to
  /// \param sigcut  Definition of signal (inverse defines background)
  /// \param presel  Cuts that will always be applied
  /// \param axes    Variables that will be cut on
  /// \param[out] cut_pos Starting cut position for each variable. One of these values will be updated
  /// \returns FOM achieved by optimized cuts
  double OptimizeOneCut(const std::string& wildcard,
                        double pot,
                        const Cut& sigcut,
                        const Cut& presel,
                        const std::vector<HistAxis>& axes,
                        std::vector<double>& cut_pos);

  /// \brief Repeatedly invoke \ref OptimizeOneCut until the FOM increase becomes small
  void OptimizeCuts(const std::string& wildcard,
                    double pot,
                    const Cut& sigcut,
                    const Cut& presel,
                    const std::vector<HistAxis>& axes,
                    std::vector<double>& cut_pos);
}
