//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/CAFAna/Core/ISyst.h"

namespace ana
{

  // Temporary pi syst

  class NuMIXSecPiSyst: public ISyst
  {
  public:

    NuMIXSecPiSyst(const std::string& name, const std::string& latexName);

    double GetSPPQ2Reweight(double Q2_GeV2, double parameter_value) const;
    double GetSPPTpiReweight(double Tpi_GeV, double parameter_value) const;

    // copied from https://github.com/GENIE-MC/Generator/blob/b336eb7b8b99957677979ab44bd0e0062f46efef/src/Framework/GHEP/GHepUtils.cxx#L22
    int NeutReactionCode(caf::SRTrueInteractionProxy *sr) const;

    // copied from https://github.com/GENIE-MC/Generator/blob/b336eb7b8b99957677979ab44bd0e0062f46efef/src/Framework/Interaction/Target.cxx#L283C14-L283C25
    bool HitNucIsSet(int pdgc) const;
    // copied from https://github.com/GENIE-MC/Generator/blob/b336eb7b8b99957677979ab44bd0e0062f46efef/src/Framework/ParticleData/PDGUtils.cxx#L346
    bool IsNucleon(int pdgc) const;
    // copied from https://github.com/GENIE-MC/Generator/blob/b336eb7b8b99957677979ab44bd0e0062f46efef/src/Framework/ParticleData/PDGUtils.cxx#L402
    bool Is2NucleonCluster(int pdgc) const;


    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:

    // pdgids
    int kPdgNeutron;
    int kPdgProton;
    int kPdgClusterNN;
    int kPdgClusterNP;
    int kPdgClusterPP;
    int kPdgPiP;
    int kPdgPiM;
    int kPdgPi0;
    int kPdgEta;
    int kPdgKP;
    int kPdgKM;
    int kPdgK0;
    int kPdgAntiK0;
    int kPdgLambda;
    int kPdgAntiLambda;
    int kPdgGamma;

  };


}
