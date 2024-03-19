//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/CAFAna/Core/ISyst.h"

namespace ana
{

  // Temporary pi syst

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

  extern const Var kNuMISPPQ2RW;

  extern const Var kNuMISPPTpiCHLinearFitReweight;
  extern const Var kNuMISPPTpiFeLinearFitReweight;
  extern const Var kNuMISPPTpiPbLinearFitReweight;

  extern const Var kNuMISPPTpiMINERvATemplateReweight;
  extern const Var kNuMISPPTpiMINERvAFittedReweight;


}
