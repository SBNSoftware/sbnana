#include "sbnana/CAFAna/Systs/EnergySysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana {

  void EnergyScaleSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    double scale = uncertainty * sigma;
    auto& det = sr->truth.det;
    bool all = (detector == EnergyScaleSystDetector::kAll);
    bool nd = (detector == EnergyScaleSystDetector::kSBND && det == caf::kSBND);
    bool ub = (detector == EnergyScaleSystDetector::kMicroBooNE && det == caf::kUNKNOWN);
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

   const EnergyScaleSyst kEnergyScaleMuon(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuon", "Correlated 2% scale systematic on E_{#mu}");
   const EnergyScaleSyst kEnergyScaleMuonSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuonSqrt", "Correlated 2% scale systematic on E_{#mu} * #sqrt{E_{#mu}}");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuonInvSqrt", "Correlated 2% scale systematic on E_{#mu} / #sqrt{E_{#mu}}");

   const EnergyScaleSyst kEnergyScaleMuonND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonND", "Uncorrelated 2% scale systematic on E_{#mu} in SBND");
   const EnergyScaleSyst kEnergyScaleMuonSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonSqrtND", "Uncorrelated 2% scale systematic on E_{#mu} * #sqrt{E_{#mu}} in SBND");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonInvSqrtND", "Uncorrelated 2% scale systematic on E_{#mu} / #sqrt{E_{#mu}} in SBND");

   const EnergyScaleSyst kEnergyScaleMuonUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonUB", "Uncorrelated 2% scale systematic on E_{#mu} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleMuonSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonSqrtUB", "Uncorrelated 2% scale systematic on E_{#mu} * #sqrt{E_{#mu}} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonInvSqrtUB", "Uncorrelated 2% scale systematic on E_{#mu} / #sqrt{E_{#mu}} in MicroBooNE");

   const EnergyScaleSyst kEnergyScaleMuonFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonFD", "Uncorrelated 2% scale systematic on E_{#mu} in ICARUS");
   const EnergyScaleSyst kEnergyScaleMuonSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonSqrtFD", "Uncorrelated 2% scale systematic on E_{#mu} * #sqrt{E_{#mu}} in ICARUS");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonInvSqrtFD", "Uncorrelated 2% scale systematic on E_{#mu} / #sqrt{E_{#mu}} in ICARUS");

   const EnergyScaleSyst kEnergyScaleHadron(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadron", "Correlated 5% scale systematic on E_{had}");
   const EnergyScaleSyst kEnergyScaleHadronSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadronSqrt", "Correlated 5% scale systematic on E_{had} * #sqrt{E_{had}}");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadronInvSqrt", "Correlated 5% scale systematic on E_{had} / #sqrt{E_{had}}");

   const EnergyScaleSyst kEnergyScaleHadronND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronND", "Uncorrelated 5% scale systematic on E_{had} in SBND");
   const EnergyScaleSyst kEnergyScaleHadronSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronSqrtND", "Uncorrelated 5% scale systematic on E_{had} * #sqrt{E_{had}} in SBND");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronInvSqrtND", "Uncorrelated 5% scale systematic on E_{had} / #sqrt{E_{had}} in SBND");

   const EnergyScaleSyst kEnergyScaleHadronUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronUB", "Uncorrelated 5% scale systematic on E_{had} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleHadronSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronSqrtUB", "Uncorrelated 5% scale systematic on E_{had} * #sqrt{E_{had}} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronInvSqrtUB", "Uncorrelated 5% scale systematic on E_{had} / #sqrt{E_{had}} in MicroBooNE");

   const EnergyScaleSyst kEnergyScaleHadronFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronFD", "Uncorrelated 5% scale systematic on E_{had} in ICARUS");
   const EnergyScaleSyst kEnergyScaleHadronSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronSqrtFD", "Uncorrelated 5% scale systematic on E_{had} * #sqrt{E_{had}} in ICARUS");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronInvSqrtFD", "Uncorrelated 5% scale systematic on E_{had} / #sqrt{E_{had}} in ICARUS");

   const EnergyScaleSyst kEnergyScaleMuonBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonBig", "Correlated 5% scale systematic on E_{#mu}");
   const EnergyScaleSyst kEnergyScaleMuonSqrtBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonSqrtBig", "Correlated 5% scale systematic on E_{#mu} * #sqrt{E_{#mu}}");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonInvSqrtBig", "Correlated 5% scale systematic on E_{#mu} / #sqrt{E_{#mu}}");

   const EnergyScaleSyst kEnergyScaleMuonNDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonNDBig", "Uncorrelated 5% scale systematic on E_{#mu} in SBND");
   const EnergyScaleSyst kEnergyScaleMuonSqrtNDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonSqrtNDBig", "Uncorrelated 5% scale systematic on E_{#mu} * #sqrt_{E_{#mu}} in SBND");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtNDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonInvSqrtNDBig", "Uncorrelated 5% scale systematic on E_{#mu} / #sqrt{E_{#mu}} in SBND");

   const EnergyScaleSyst kEnergyScaleMuonUBBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonUBBig", "Uncorrelated 5% scale systematic on E_{#mu} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleMuonSqrtUBBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonSqrtUBBig", "Uncorrelated 5% scale systematic on E_{#mu} * #sqrt{E_{#mu}} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtUBBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonInvSqrtUBBig", "Uncorrelated 5% scale systematic on E_{#mu} / #sqrt{E_{#mu}} in MicroBooNE");

   const EnergyScaleSyst kEnergyScaleMuonFDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonFDBig", "Uncorrelated 5% scale systematic on E_{#mu} in ICARUS");
   const EnergyScaleSyst kEnergyScaleMuonSqrtFDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonSqrtFDBig", "Uncorrelated 5% scale systematic on E_{#mu} * #sqrt{E_{#mu}} in ICARUS");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtFDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonInvSqrtFDBig", "Uncorrelated 5% scale systematic on E_{#mu} / #sqrt{E_{#mu}} in ICARUS");

   const EnergyScaleSyst kEnergyScaleHadronBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronBig", "Correlated 10% scale systematic on E_{had}");
   const EnergyScaleSyst kEnergyScaleHadronSqrtBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronSqrtBig", "Correlated 10% scale systematic on E_{had} * #sqrt{E_{had}}");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronInvSqrtBig", "Correlated 10% scale systematic on E_{had} / #sqrt{E_{had}}");

   const EnergyScaleSyst kEnergyScaleHadronNDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronNDBig", "Uncorrelated 10% scale systematic on E_{had} in SBND");
   const EnergyScaleSyst kEnergyScaleHadronSqrtNDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronSqrtNDBig", "Uncorrelated 10% scale systematic on E_{had} * #sqrt{E_{had}} in SBND");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtNDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronInvSqrtNDBig", "Uncorrelated 10% scale systematic on E_{had} / #sqrt{E_{had}} in SBND");

   const EnergyScaleSyst kEnergyScaleHadronUBBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronUBBig", "Uncorrelated 10% scale systematic on E_{had} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleHadronSqrtUBBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronSqrtUBBig", "Uncorrelated 10% scale systematic on E_{had} * #sqrt{E_had}} in MicroBooNE");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtUBBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronInvSqrtUBBig", "Uncorrelated 10% scale systematic on E_{had} / #sqrt{E_{had}} in MicroBooNE");

   const EnergyScaleSyst kEnergyScaleHadronFDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronFDBig", "Uncorrelated 10% scale systematic on E_{had} in ICARUS");
   const EnergyScaleSyst kEnergyScaleHadronSqrtFDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronSqrtFDBig", "Uncorrelated 10% scale systematic on E_{had} * #sqrt{E_{had}} in ICARUS");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtFDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronInvSqrtFDBig", "Uncorrelated 10% scale systematic on E_{had} / #sqrt{E_{had}} in ICARUS");

  std::vector<const EnergyScaleSyst*> GetEnergySysts() {
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

  std::vector<const EnergyScaleSyst*> GetBigEnergySysts() {
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
