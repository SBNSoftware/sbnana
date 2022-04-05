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

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part, EnergyScaleSystDetector detector = EnergyScaleSystDetector::kAll>
  class EnergyScaleSyst: public ISyst
  {
  public:
    EnergyScaleSyst(double _uncertainty, std::string name) : 
      ISyst(name, name), uncertainty(_uncertainty) {}
    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override
    {
      double scale = uncertainty * sigma;
      auto baseline_cut = [](double baseline){ 
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
    double uncertainty;
  };

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part>
  using EnergyScaleSystCorr = EnergyScaleSyst<term, part, EnergyScaleSystDetector::kAll>;

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part>
  using EnergyScaleSystUncorrND = EnergyScaleSyst<term, part, EnergyScaleSystDetector::kSBND>;

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part>
  using EnergyScaleSystUncorrUB = EnergyScaleSyst<term, part, EnergyScaleSystDetector::kMicroBooNE>;

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part>
  using EnergyScaleSystUncorrFD = EnergyScaleSyst<term, part, EnergyScaleSystDetector::kICARUS>;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuon;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonSqrt;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonInvSqrt;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtND;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtUB;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtFD;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadron;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronSqrt;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronInvSqrt;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtND;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtUB;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtFD;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonSqrtBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonInvSqrtBig;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtNDBig;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtUBBig;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtFDBig;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronSqrtBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronInvSqrtBig;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtNDBig;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtUBBig;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtFDBig;

  // Vector of energy scale systematics
  struct EnergySystVector: public std::vector<const ISyst*>
  {

  };

EnergySystVector GetEnergySysts();
EnergySystVector GetBigEnergySysts();

}
