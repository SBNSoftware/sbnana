//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"

#include <cassert>

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
    EnergyScaleSyst(EnergyScaleSystTerm _term, EnergyScaleSystParticle _part, EnergyScaleSystDetector _detector, double _uncertainty, std::string name) : 
      ISyst(name, name), term(_term), part(_part), detector(_detector), uncertainty(_uncertainty) {}
    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override
    {
      double scale = uncertainty * sigma;
      auto baseline_cut = [detector = detector](double baseline){ 
        auto all = (detector == EnergyScaleSystDetector::kAll);
        auto nd = (detector == EnergyScaleSystDetector::kSBND && baseline < 120);
        auto ub = (detector == EnergyScaleSystDetector::kMicroBooNE && baseline > 120 && baseline < 500);
        auto fd = (detector == EnergyScaleSystDetector::kICARUS && baseline > 500);
        return all || nd || ub || fd; 
      };
      if(sr->truth.iscc && abs(sr->truth.pdg) == 14 && !isnan(sr->fake_reco.nuE) && baseline_cut(sr->truth.baseline)) {
        double particle_energy = 0.0;
        switch(part) {
        case EnergyScaleSystParticle::kMuon:
          particle_energy += sr->fake_reco.lepton.ke;
          break;
        case EnergyScaleSystParticle::kHadron:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            particle_energy += sr->fake_reco.hadrons[i].ke;
          }
          break;
        case EnergyScaleSystParticle::kNeutron:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            if(sr->fake_reco.hadrons[i].pid == 2112) {
              particle_energy += sr->fake_reco.hadrons[i].ke;
            }
          }
          break;
        case EnergyScaleSystParticle::kEM:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            if(sr->fake_reco.hadrons[i].pid == 111) {
              particle_energy += sr->fake_reco.hadrons[i].ke;
            }
          }
          break;
        case EnergyScaleSystParticle::kChargedHadron:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            auto pid = sr->fake_reco.hadrons[i].pid;
            if(pid == 2212 || abs(pid) == 211) {
              particle_energy += sr->fake_reco.hadrons[i].ke;
            }
          }
          break;
        }
        switch(term) {
        case EnergyScaleSystTerm::kConstant:
          break;
        case EnergyScaleSystTerm::kSqrt:
          scale *= std::sqrt(particle_energy);
          break;
        case EnergyScaleSystTerm::kInverseSqrt:
          scale /= std::sqrt(particle_energy + 0.1);
          break;
        }
        sr->fake_reco.nuE += particle_energy * scale;
      }
    }
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

  // Vector of energy scale systematics
  struct EnergySystVector: public std::vector<const ISyst*>
  {

  };

EnergySystVector GetEnergySysts();
EnergySystVector GetBigEnergySysts();

}
