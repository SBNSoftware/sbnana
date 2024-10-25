//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include <vector>

namespace ana
{

  enum class EnergyScaleSystTerm {
    kConstant,
    kSqrt,
    kInverseSqrt
  };

  enum class EnergyScaleSystParticle {
    kMuon,
    kHadron,
    kNeutron,
    kEM,
    kChargedHadron
  };

  enum class EnergyScaleSystDetector {
    kSBND,
    kMicroBooNE,
    kICARUS,
    kAll
  };

  class EnergyScaleSyst: public ISyst
  {
  public:
    EnergyScaleSyst(EnergyScaleSystTerm _term, EnergyScaleSystParticle _part, EnergyScaleSystDetector _detector, double _uncertainty, const std::string& name, const std::string& latexName) : 
      ISyst(name, latexName), term(_term), part(_part), detector(_detector), uncertainty(_uncertainty) {}

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:
    EnergyScaleSystTerm term;
    EnergyScaleSystParticle part;
    EnergyScaleSystDetector detector;
    double uncertainty;
  };

  class RecoEnergyScaleSyst: public ISyst
  {
  public:
    RecoEnergyScaleSyst(EnergyScaleSystTerm _term, EnergyScaleSystParticle _part, EnergyScaleSystDetector _detector, double _uncertainty, const std::string& name, const std::string& latexName) : 
      ISyst(name, latexName), term(_term), part(_part), detector(_detector), uncertainty(_uncertainty) {}

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:
    EnergyScaleSystTerm term;
    EnergyScaleSystParticle part;
    EnergyScaleSystDetector detector;
    double uncertainty;
  };

  extern const EnergyScaleSyst kEnergyScaleMuon;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrt;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrt;

  extern const EnergyScaleSyst kEnergyScaleMuonND;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrtND;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrtND;

  extern const EnergyScaleSyst kEnergyScaleMuonUB;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrtUB;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrtUB;

  extern const EnergyScaleSyst kEnergyScaleMuonFD;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrtFD;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrtFD;

  extern const EnergyScaleSyst kEnergyScaleHadron;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrt;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrt;

  extern const EnergyScaleSyst kEnergyScaleHadronND;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrtND;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrtND;

  extern const EnergyScaleSyst kEnergyScaleHadronUB;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrtUB;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrtUB;

  extern const EnergyScaleSyst kEnergyScaleHadronFD;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrtFD;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrtFD;

  extern const EnergyScaleSyst kEnergyScaleMuonBig;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrtBig;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrtBig;

  extern const EnergyScaleSyst kEnergyScaleMuonNDBig;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrtNDBig;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrtNDBig;

  extern const EnergyScaleSyst kEnergyScaleMuonUBBig;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrtUBBig;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrtUBBig;

  extern const EnergyScaleSyst kEnergyScaleMuonFDBig;
  extern const EnergyScaleSyst kEnergyScaleMuonSqrtFDBig;
  extern const EnergyScaleSyst kEnergyScaleMuonInvSqrtFDBig;

  extern const EnergyScaleSyst kEnergyScaleHadronBig;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrtBig;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrtBig;

  extern const EnergyScaleSyst kEnergyScaleHadronNDBig;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrtNDBig;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrtNDBig;

  extern const EnergyScaleSyst kEnergyScaleHadronUBBig;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrtUBBig;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrtUBBig;

  extern const EnergyScaleSyst kEnergyScaleHadronFDBig;
  extern const EnergyScaleSyst kEnergyScaleHadronSqrtFDBig;
  extern const EnergyScaleSyst kEnergyScaleHadronInvSqrtFDBig;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuon;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrt;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrt;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonND;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrtND;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrtND;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonUB;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrtUB;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrtUB;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonFD;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrtFD;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrtFD;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadron;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrt;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrt;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronND;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrtND;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrtND;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronUB;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrtUB;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrtUB;

  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronFD;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrtFD;
  extern const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrtFD;

  std::vector<const ISyst*> GetEnergySysts();
  std::vector<const ISyst*> GetBigEnergySysts();

  std::vector<const ISyst*> GetRecoEnergySysts();

}
