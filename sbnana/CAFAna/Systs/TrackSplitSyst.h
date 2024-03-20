// Bruce Howard based on CalorimetrySysts -- 2024
//
// SEE THE .CXX FILE FOR AN IMPORTANT NOTE ON WHAT'S MODIFIED OR NOT.

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "TFile.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSpline.h"

#include <vector>

namespace ana
{

  struct TrkChi2Results { ///< determined particle ID
    Float_t chi2_kaon, chi2_muon, chi2_pion, chi2_proton, pida;
    Int_t pid_ndof;
  };

  struct TrkMomentumResults { ///< range-based momentum results
    Float_t p_muon, p_pion, p_proton;
  };

  struct PreservedInitialTrkResults { ///< stuff from the unchanged initial track that I want to copy to the secondary track
    Float_t dir_end_x, dir_end_y, dir_end_z, end_x, end_y, end_z;
    unsigned producer;
  };

  class TrackSplitSyst: public ISyst
  {
  public:

    TrackSplitSyst(const std::string& name, const std::string& latexName, const bool debugPrint=false, const std::string& fileName="");

    TrkChi2Results CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo) const;
    TrkChi2Results CalculateChi2(const caf::SRTrackCalo& calo) const;
    TrkMomentumResults CalculateMomenta(const float length) const;

    caf::SRCaloPoint FillCaloPointFrom( const caf::Proxy<caf::SRCaloPoint>& inCaloPt ) const;
    void FillPtrPFP( caf::SRPFP& ret, const caf::Proxy<caf::SRPFP>& inPfp ) const;

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:
    // Random engine
    TRandom3 *tRand;

    // split probabilities
    TH2D *splitProb[2][2]; // In ICARUS, idx1: 0 = East and 1 = West, idx2: 0 = cathode, 1 = z=0 point

    // Externals
    TSpline3 KEvsR_spline3;

    TProfile *dedx_range_pro;   ///< proton template
    TProfile *dedx_range_ka;    ///< kaon template
    TProfile *dedx_range_pi;    ///< pion template
    TProfile *dedx_range_mu;    ///< muon template

    // Debug printouts?
    bool fDebug;
  };

  extern const TrackSplitSyst kTrackSplittingSyst;
  extern const TrackSplitSyst kTrackSplittingSystDebug;


  class TrackSplitSystCheck: public TrackSplitSyst
  {
  public:

    TrackSplitSystCheck(const std::string& name, const std::string& latexName, const bool debugPrint=false, const std::string& fileName="")
    : TrackSplitSyst( name, latexName, debugPrint, fileName ) { }

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;
  };

  extern const TrackSplitSystCheck kTrackSplittingSystCheck;

}
