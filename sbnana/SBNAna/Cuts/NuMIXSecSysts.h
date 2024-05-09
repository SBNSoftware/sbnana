//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/CAFAna/Core/ISyst.h"

#include "TH1.h"

namespace ana
{

  //---------------------------------------------------------------------
  // Single-pion production CV correction and systematics

  bool IsSPP(const caf::SRTrueInteractionProxy *sr);
  extern const Var kNuMITrueIsSPP;

  // Q2 template RW
  double GetSPPQ2Reweight(double Q2_GeV2);
  // Tpi RW
  // - CH, linear fit
  double GetSPPTpiCHLinearFitReweight(double Tpi_GeV);
  // - Fe, linear fit
  double GetSPPTpiFeLinearFitReweight(double Tpi_GeV);
  // - Pb, linear fit
  double GetSPPTpiPbLinearFitReweight(double Tpi_GeV);
  // - MINERvA, tempalte
  double GetSPPTpiMINERvATemplateReweight(double Tpi_GeV);
  // - MINERvA, fitted (tpi<225MeV), linear interplolation (225<=tpi<237.5MeV), templtae (tpi>=237.5MEV)
  double GetSPPTpiMINERvAFittedReweight(double Tpi_GeV);

  class NuMIXSecPiSyst: public ISyst
  {
  public:

    NuMIXSecPiSyst(const std::string& name, const std::string& latexName);

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:

  };

  // CV correction weight; kNuMISPPQ2RW * kNuMISPPTpiMINERvAFittedReweight
  extern const TruthVar kTruth_NuMISPPCVCorrection;
  extern const Var kNuMISPPCVCorrection;

  // Separate reweight values for study
  // - Q2-template RW
  extern const Var kNuMISPPQ2RW;
  // - Tpi linear-fitted RW using MINERvA CC1pi; data is available for 35<tpi<350 MeV,
  //   and the RW is extrapolated down to 0. Above 350 MeV, the value at 350 MeV is used
  extern const Var kNuMISPPTpiCHLinearFitReweight;
  extern const Var kNuMISPPTpiFeLinearFitReweight;
  extern const Var kNuMISPPTpiPbLinearFitReweight;
  // - MINERvA untracked pion analysis observed a tpi suppresion for low tpi from data
  //   - binned value
  extern const Var kNuMISPPTpiMINERvATemplateReweight;
  //   - binned value fitted with analytic function
  extern const Var kNuMISPPTpiMINERvAFittedReweight;

  //---------------------------------------------------------------------
  // Split-track reweighting

  class NuMIXSecSplitTrackReweight{

    public:

      NuMIXSecSplitTrackReweight();
      ~NuMIXSecSplitTrackReweight();
      static NuMIXSecSplitTrackReweight& Instance();

      TH1* fRWCathode[2][2]; // [cryo; 0/1][IsSplit?]

      int CathodeSplitType(const caf::Proxy<caf::SRTrack>& trk) const;
      double GetCathodeRW(const caf::Proxy<caf::SRTrack>& trk) const;

  };

  // CV correction
  extern const Var kNuMISplitTrackCVCorrection;

/*
  class NuMIXSecSplitTrackSyst: public ISyst
  {
  public:

    NuMIXSecSplitTrackSyst(const std::string& name, const std::string& latexName);

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:

  };
*/

}
