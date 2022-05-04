#include "sbnana/CAFAna/Systs/EnergySysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana {

  void EnergyScaleSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    double scale = uncertainty * sigma;
    auto& det = sr->truth.det;
    bool all = (detector == EnergyScaleSystDetector::kAll);
    bool nd = (detector == EnergyScaleSystDetector::kSBND && det == caf::kSBND);
    bool ub = (detector == EnergyScaleSystDetector::kMicroBooNE && det != caf::kSBND && det != caf::kICARUS);
    bool fd = (detector == EnergyScaleSystDetector::kICARUS && det == caf::kICARUS);
    bool detector_cut = all || nd || ub || fd;
    if(!sr->truth.iscc || abs(sr->truth.pdg) != 14 || isnan(sr->fake_reco.nuE) || !detector_cut) 
      return ;
    double particle_energy = 0.0;
    switch(part) {
    case EnergyScaleSystParticle::kMuon:
      particle_energy += sr->fake_reco.lepton.ke;
      break;
    case EnergyScaleSystParticle::kHadron:
      for(const auto& hadron: sr->fake_reco.hadrons) {
        particle_energy += hadron.ke;
      }
      break;
    case EnergyScaleSystParticle::kNeutron:
      for(const auto& hadron: sr->fake_reco.hadrons) {
        if(hadron.pid == 2112) {
          particle_energy += hadron.ke;
        }
      }
      break;
    case EnergyScaleSystParticle::kEM:
      for(const auto& hadron: sr->fake_reco.hadrons) {
        if(hadron.pid == 111) {
          particle_energy += hadron.ke;
        }
      }
      break;
   case EnergyScaleSystParticle::kChargedHadron:
      for(const auto& hadron: sr->fake_reco.hadrons) {
        auto pid = hadron.pid;
        if(pid == 2212 || abs(pid) == 211) {
          particle_energy += hadron.ke;
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

   const EnergyScaleSyst kEnergyScaleMuon(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuon", "Correlated linear E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuonSqrt", "Correlated sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuonInvSqrt", "Correlated inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleMuonND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonND", "Uncorrelated SBND linear E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonSqrtND", "Uncorrelated SBND sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonInvSqrtND", "Uncorrelated SBND inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleMuonUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonUB", "Uncorrelated MicroBooNE linear E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonSqrtUB", "Uncorrelated MicroBooNE sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonInvSqrtUB", "Uncorrelated MicroBooNE inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleMuonFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonFD", "Uncorrelated ICARUS linear E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonSqrtFD", "Uncorrelated ICARUS sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonInvSqrtFD", "Uncorrelated ICARUS inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleHadron(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadron", "Correlated linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadronSqrt", "Correlated sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadronInvSqrt", "Correlated inv sqrt E_{had} scale");

   const EnergyScaleSyst kEnergyScaleHadronND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronND", "Uncorrelated SBND linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronSqrtND", "Uncorrelated SBND sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronInvSqrtND", "Uncorrelated SBND inv sqrt E_{had} scale");

   const EnergyScaleSyst kEnergyScaleHadronUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronUB", "Uncorrelated MicroBooNE linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronSqrtUB", "Uncorrelated MicroBooNE sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronInvSqrtUB", "Uncorrelated MicroBooNE inv sqrt E_{had} scale");

   const EnergyScaleSyst kEnergyScaleHadronFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronFD", "Uncorrelated ICARUS linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronSqrtFD", "Uncorrelated ICARUS sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronInvSqrtFD", "Uncorrelated ICARUS inv sqrt E_{had} scale");

   const EnergyScaleSyst kEnergyScaleMuonBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonBig", "Correlated linear E_{#mu} scale ");
   const EnergyScaleSyst kEnergyScaleMuonSqrtBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonSqrtBig", "Correlated sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonInvSqrtBig", "Correlated inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleMuonNDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonNDBig", "Uncorrelated SBND linear E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonSqrtNDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonSqrtNDBig", "Uncorrelated SBND sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtNDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonInvSqrtNDBig", "Uncorrelated SBND inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleMuonUBBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonUBBig", "Uncorrelated MicroBooNE linear E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonSqrtUBBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonSqrtUBBig", "Uncorrelated MicroBooNE sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtUBBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonInvSqrtUBBig", "Uncorrelated MicroBooNE inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleMuonFDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonFDBig", "Uncorrelated ICARUS linear E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonSqrtFDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonSqrtFDBig", "Uncorrelated ICARUS sqrt E_{#mu} scale");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtFDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonInvSqrtFDBig", "Uncorrelated ICARUS inv sqrt E_{#mu} scale");

   const EnergyScaleSyst kEnergyScaleHadronBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronBig", "Correlated linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrtBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronSqrtBig", "Correlated sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronInvSqrtBig", "Correlated inv sqrt E_{had} scale");

   const EnergyScaleSyst kEnergyScaleHadronNDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronNDBig", "Uncorrelated SBND linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrtNDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronSqrtNDBig", "Uncorrelated SBND sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtNDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronInvSqrtNDBig", "Uncorrelated SBND inv sqrt E_{had} scale");

   const EnergyScaleSyst kEnergyScaleHadronUBBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronUBBig", "Uncorrelated MicroBooNE linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrtUBBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronSqrtUBBig", "Uncorrelated MicroBooNE sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtUBBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronInvSqrtUBBig", "Uncorrelated MicroBooNE inv sqrt E_{had} scale");

   const EnergyScaleSyst kEnergyScaleHadronFDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronFDBig", "Uncorrelated ICARUS linear E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronSqrtFDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronSqrtFDBig", "Uncorrelated ICARUS sqrt E_{had} scale");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtFDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronInvSqrtFDBig", "Uncorrelated ICARUS inv sqrt E_{had} scale");

  std::vector<const ISyst*> GetEnergySysts() {
    // MicroBooNE, ChargedHadron, Neutron, and EM systs
    // not being used right now
    return {&kEnergyScaleMuon,
            &kEnergyScaleMuonSqrt,
            &kEnergyScaleMuonInvSqrt,
            &kEnergyScaleMuonND,
            &kEnergyScaleMuonSqrtND,
            &kEnergyScaleMuonInvSqrtND,
            //&kEnergyScaleMuonUB,
            //&kEnergyScaleMuonSqrtUB,
            //&kEnergyScaleMuonInvSqrtUB,
            &kEnergyScaleMuonFD,
            &kEnergyScaleMuonSqrtFD,
            &kEnergyScaleMuonInvSqrtFD,
            &kEnergyScaleHadron,
            //&kEnergyScaleChargedHadron,
            //&kEnergyScaleEM,
            //&kEnergyScaleNeutron,
            &kEnergyScaleHadronSqrt,
            //&kEnergyScaleEMSqrt,
            //&kEnergyScaleNeutronSqrt,
            &kEnergyScaleHadronInvSqrt,
            //&kEnergyScaleEMInvSqrt,
            //&kEnergyScaleNeutronInvSqrt,
            &kEnergyScaleHadronND,
            //&kEnergyScaleChargedHadronND,
            //&kEnergyScaleEMND,
            //&kEnergyScaleNeutronND,
            &kEnergyScaleHadronSqrtND,
            //&kEnergyScaleEMSqrtND,
            //&kEnergyScaleNeutronSqrtND,
            &kEnergyScaleHadronInvSqrtND,
            //&kEnergyScaleEMInvSqrtND,
            //&kEnergyScaleNeutronInvSqrtND,
            //&kEnergyScaleHadronUB,
            //&kEnergyScaleChargedHadronUB,
            //&kEnergyScaleEMUB,
            //&kEnergyScaleNeutronUB,
            //&kEnergyScaleHadronSqrtUB,
            //&kEnergyScaleEMSqrtUB,
            //&kEnergyScaleNeutronSqrtUB,
            //&kEnergyScaleHadronInvSqrtUB,
            //&kEnergyScaleEMInvSqrtUB,
            //&kEnergyScaleNeutronInvSqrtUB,
            &kEnergyScaleHadronFD,
            //&kEnergyScaleChargedHadronFD,
            //&kEnergyScaleEMFD,
            //&kEnergyScaleNeutronFD,
            &kEnergyScaleHadronSqrtFD,
            //&kEnergyScaleEMSqrtFD,
            //&kEnergyScaleNeutronSqrtFD,
            //&kEnergyScaleEMInvSqrtFD,
            //&kEnergyScaleNeutronInvSqrtFD,
            &kEnergyScaleHadronInvSqrtFD};
  }

  std::vector<const ISyst*> GetBigEnergySysts() {
    return {&kEnergyScaleMuonBig,
            &kEnergyScaleMuonSqrtBig,
            &kEnergyScaleMuonInvSqrtBig,
            &kEnergyScaleMuonNDBig,
            &kEnergyScaleMuonSqrtNDBig,
            &kEnergyScaleMuonInvSqrtNDBig,
            //&kEnergyScaleMuonUBBig,
            //&kEnergyScaleMuonSqrtUBBig,
            //&kEnergyScaleMuonInvSqrtUBBig,
            &kEnergyScaleMuonFDBig,
            &kEnergyScaleMuonSqrtFDBig,
            &kEnergyScaleMuonInvSqrtFDBig,
            &kEnergyScaleHadronBig,
            &kEnergyScaleHadronSqrtBig,
            &kEnergyScaleHadronInvSqrtBig,
            &kEnergyScaleHadronNDBig,
            &kEnergyScaleHadronSqrtNDBig,
            &kEnergyScaleHadronInvSqrtNDBig,
            //&kEnergyScaleHadronUBBig,
            //&kEnergyScaleHadronSqrtUBBig,
            //&kEnergyScaleHadronInvSqrtUBBig,
            &kEnergyScaleHadronFDBig,
            &kEnergyScaleHadronSqrtFDBig,
            &kEnergyScaleHadronInvSqrtFDBig};
  }
} // namespace ana
