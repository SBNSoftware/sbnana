//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"

#include <cassert>

namespace ana
{

  enum class EnergyScaleSystTerm {
    Constant,
    Sqrt,
    InverseSqrt
  };

  enum class EnergyScaleSystParticle {
    Muon,
    Hadron,
    Neutron,
    Pi0,
    ChargedHadron
  };

  enum class EnergyScaleSystDetector {
    SBND,
    MicroBooNE,
    ICARUS,
    All
  };

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part, EnergyScaleSystDetector detector = EnergyScaleSystDetector::All>
  class EnergyScaleSyst: public ISyst
  {
  public:
    EnergyScaleSyst(double _uncertainty, std::string name) : 
      ISyst(name, name), uncertainty(_uncertainty) {}
    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override
    {
      double scale = uncertainty * sigma;
      auto baseline_cut = [](double baseline){ 
        auto all = (detector == EnergyScaleSystDetector::All);
        auto nd = (detector == EnergyScaleSystDetector::SBND && baseline < 120);
        auto ub = (detector == EnergyScaleSystDetector::MicroBooNE && baseline > 800);
        auto fd = (detector == EnergyScaleSystDetector::ICARUS && baseline > 500);
        return all || nd || ub || fd; 
      };
      if(sr->truth.iscc && abs(sr->truth.pdg) == 14 && !isnan(sr->fake_reco.nuE) && baseline_cut(sr->truth.baseline)) {
        double particle_energy = 0.0;
        switch(part) {
        case EnergyScaleSystParticle::Muon:
          particle_energy += sr->fake_reco.lepton.ke;
          break;
        case EnergyScaleSystParticle::Hadron:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            particle_energy += sr->fake_reco.hadrons[i].ke;
          }
          break;
        case EnergyScaleSystParticle::Neutron:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            if(sr->fake_reco.hadrons[i].pid == 2112) {
              particle_energy += sr->fake_reco.hadrons[i].ke;
            }
          }
          break;
        case EnergyScaleSystParticle::Pi0:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            if(sr->fake_reco.hadrons[i].pid == 111) {
              particle_energy += sr->fake_reco.hadrons[i].ke;
            }
          }
          break;
        case EnergyScaleSystParticle::ChargedHadron:
          for (size_t i = 0; i < sr->fake_reco.hadrons.size(); ++i) {
            auto pid = sr->fake_reco.hadrons[i].pid;
            if(pid == 2212 || abs(pid) == 211) {
              particle_energy += sr->fake_reco.hadrons[i].ke;
            }
          }
          break;
        }
        switch(term) {
        case EnergyScaleSystTerm::Constant:
          break;
        case EnergyScaleSystTerm::Sqrt:
          scale *= std::sqrt(particle_energy);
          break;
        case EnergyScaleSystTerm::InverseSqrt:
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
  using EnergyScaleSystCorr = EnergyScaleSyst<term, part, EnergyScaleSystDetector::All>;

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part>
  using EnergyScaleSystUncorrND = EnergyScaleSyst<term, part, EnergyScaleSystDetector::SBND>;

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part>
  using EnergyScaleSystUncorrUB = EnergyScaleSyst<term, part, EnergyScaleSystDetector::MicroBooNE>;

  template<EnergyScaleSystTerm term, EnergyScaleSystParticle part>
  using EnergyScaleSystUncorrFD = EnergyScaleSyst<term, part, EnergyScaleSystDetector::ICARUS>;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuon;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonSqrt;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonInvSqrt;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtND;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtUB;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtFD;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadron;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronSqrt;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronInvSqrt;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtND;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtND;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtUB;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtUB;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtFD;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtFD;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonSqrtBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonInvSqrtBig;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtNDBig;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtUBBig;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtFDBig;

  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronSqrtBig;
  extern const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronInvSqrtBig;

  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtNDBig;
  extern const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtNDBig;

  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtUBBig;
  extern const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtUBBig;

  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtFDBig;
  extern const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtFDBig;

  // Vector of energy scale systematics
  struct EnergySystVector: public std::vector<const ISyst*>
  {

  };

EnergySystVector GetEnergySysts();
EnergySystVector GetBigEnergySysts();

}
