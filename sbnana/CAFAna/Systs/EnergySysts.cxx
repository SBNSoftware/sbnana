#include "sbnana/CAFAna/Systs/EnergySysts.h"

namespace ana {

  void EnergyScaleSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
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

   const EnergyScaleSyst kEnergyScaleMuon(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuon");
   const EnergyScaleSyst kEnergyScaleMuonSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuonSqrt");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "EnergyScaleMuonInvSqrt");

   const EnergyScaleSyst kEnergyScaleMuonND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonND");
   const EnergyScaleSyst kEnergyScaleMuonSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonSqrtND");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "EnergyScaleMuonInvSqrtND");

   const EnergyScaleSyst kEnergyScaleMuonUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonUB");
   const EnergyScaleSyst kEnergyScaleMuonSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonSqrtUB");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "EnergyScaleMuonInvSqrtUB");

   const EnergyScaleSyst kEnergyScaleMuonFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonFD");
   const EnergyScaleSyst kEnergyScaleMuonSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonSqrtFD");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "EnergyScaleMuonInvSqrtFD");

   const EnergyScaleSyst kEnergyScaleHadron(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadron");
   const EnergyScaleSyst kEnergyScaleHadronSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadronSqrt");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleHadronInvSqrt");

   const EnergyScaleSyst kEnergyScaleHadronND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronND");
   const EnergyScaleSyst kEnergyScaleHadronSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronSqrtND");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleHadronInvSqrtND");

   const EnergyScaleSyst kEnergyScaleHadronUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronUB");
   const EnergyScaleSyst kEnergyScaleHadronSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronSqrtUB");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleHadronInvSqrtUB");

   const EnergyScaleSyst kEnergyScaleHadronFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronFD");
   const EnergyScaleSyst kEnergyScaleHadronSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronSqrtFD");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleHadronInvSqrtFD");

   const EnergyScaleSyst kEnergyScaleMuonBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonBig");
   const EnergyScaleSyst kEnergyScaleMuonSqrtBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonSqrtBig");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.05, "EnergyScaleMuonInvSqrtBig");

   const EnergyScaleSyst kEnergyScaleMuonNDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonNDBig");
   const EnergyScaleSyst kEnergyScaleMuonSqrtNDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonSqrtNDBig");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtNDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.05, "EnergyScaleMuonInvSqrtNDBig");

   const EnergyScaleSyst kEnergyScaleMuonUBBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonUBBig");
   const EnergyScaleSyst kEnergyScaleMuonSqrtUBBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonSqrtUBBig");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtUBBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.05, "EnergyScaleMuonInvSqrtUBBig");

   const EnergyScaleSyst kEnergyScaleMuonFDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonFDBig");
   const EnergyScaleSyst kEnergyScaleMuonSqrtFDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonSqrtFDBig");
   const EnergyScaleSyst kEnergyScaleMuonInvSqrtFDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.05, "EnergyScaleMuonInvSqrtFDBig");

   const EnergyScaleSyst kEnergyScaleHadronBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronBig");
   const EnergyScaleSyst kEnergyScaleHadronSqrtBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronSqrtBig");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.10, "EnergyScaleHadronInvSqrtBig");

   const EnergyScaleSyst kEnergyScaleHadronNDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronNDBig");
   const EnergyScaleSyst kEnergyScaleHadronSqrtNDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronSqrtNDBig");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtNDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.10, "EnergyScaleHadronInvSqrtNDBig");

   const EnergyScaleSyst kEnergyScaleHadronUBBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronUBBig");
   const EnergyScaleSyst kEnergyScaleHadronSqrtUBBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronSqrtUBBig");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtUBBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.10, "EnergyScaleHadronInvSqrtUBBig");

   const EnergyScaleSyst kEnergyScaleHadronFDBig(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronFDBig");
   const EnergyScaleSyst kEnergyScaleHadronSqrtFDBig(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronSqrtFDBig");
   const EnergyScaleSyst kEnergyScaleHadronInvSqrtFDBig(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.10, "EnergyScaleHadronInvSqrtFDBig");

  EnergySystVector GetEnergySysts() {
    // MicroBooNE, ChargedHadron, Neutron, and EM systs
    // not being used right now
    EnergySystVector vec;
    vec.push_back(&kEnergyScaleMuon);
    vec.push_back(&kEnergyScaleMuonSqrt);
    vec.push_back(&kEnergyScaleMuonInvSqrt);
    vec.push_back(&kEnergyScaleMuonND);
    vec.push_back(&kEnergyScaleMuonSqrtND);
    vec.push_back(&kEnergyScaleMuonInvSqrtND);
    //vec.push_back(&kEnergyScaleMuonUB);
    //vec.push_back(&kEnergyScaleMuonSqrtUB);
    //vec.push_back(&kEnergyScaleMuonInvSqrtUB);
    vec.push_back(&kEnergyScaleMuonFD);
    vec.push_back(&kEnergyScaleMuonSqrtFD);
    vec.push_back(&kEnergyScaleMuonInvSqrtFD);
    vec.push_back(&kEnergyScaleHadron);
    //vec.push_back(&kEnergyScaleChargedHadron);
    //vec.push_back(&kEnergyScaleEM);
    //vec.push_back(&kEnergyScaleNeutron);
    vec.push_back(&kEnergyScaleHadronSqrt);
    //vec.push_back(&kEnergyScaleEMSqrt);
    //vec.push_back(&kEnergyScaleNeutronSqrt);
    vec.push_back(&kEnergyScaleHadronInvSqrt);
    //vec.push_back(&kEnergyScaleEMInvSqrt);
    //vec.push_back(&kEnergyScaleNeutronInvSqrt);
    vec.push_back(&kEnergyScaleHadronND);
    //vec.push_back(&kEnergyScaleChargedHadronND);
    //vec.push_back(&kEnergyScaleEMND);
    //vec.push_back(&kEnergyScaleNeutronND);
    vec.push_back(&kEnergyScaleHadronSqrtND);
    //vec.push_back(&kEnergyScaleEMSqrtND);
    //vec.push_back(&kEnergyScaleNeutronSqrtND);
    vec.push_back(&kEnergyScaleHadronInvSqrtND);
    //vec.push_back(&kEnergyScaleEMInvSqrtND);
    //vec.push_back(&kEnergyScaleNeutronInvSqrtND);
    //vec.push_back(&kEnergyScaleHadronUB);
    //vec.push_back(&kEnergyScaleChargedHadronUB);
    //vec.push_back(&kEnergyScaleEMUB);
    //vec.push_back(&kEnergyScaleNeutronUB);
    //vec.push_back(&kEnergyScaleHadronSqrtUB);
    //vec.push_back(&kEnergyScaleEMSqrtUB);
    //vec.push_back(&kEnergyScaleNeutronSqrtUB);
    //vec.push_back(&kEnergyScaleHadronInvSqrtUB);
    //vec.push_back(&kEnergyScaleEMInvSqrtUB);
    //vec.push_back(&kEnergyScaleNeutronInvSqrtUB);
    vec.push_back(&kEnergyScaleHadronFD);
    //vec.push_back(&kEnergyScaleChargedHadronFD);
    //vec.push_back(&kEnergyScaleEMFD);
    //vec.push_back(&kEnergyScaleNeutronFD);
    vec.push_back(&kEnergyScaleHadronSqrtFD);
    //vec.push_back(&kEnergyScaleEMSqrtFD);
    //vec.push_back(&kEnergyScaleNeutronSqrtFD);
    vec.push_back(&kEnergyScaleHadronInvSqrtFD);
    //vec.push_back(&kEnergyScaleEMInvSqrtFD);
    //vec.push_back(&kEnergyScaleNeutronInvSqrtFD);

    return vec;
  }

  EnergySystVector GetBigEnergySysts() {

    EnergySystVector vec;
    vec.push_back(&kEnergyScaleMuonBig);
    vec.push_back(&kEnergyScaleMuonSqrtBig);
    vec.push_back(&kEnergyScaleMuonInvSqrtBig);
    vec.push_back(&kEnergyScaleMuonNDBig);
    vec.push_back(&kEnergyScaleMuonSqrtNDBig);
    vec.push_back(&kEnergyScaleMuonInvSqrtNDBig);
    //vec.push_back(&kEnergyScaleMuonUBBig);
    //vec.push_back(&kEnergyScaleMuonSqrtUBBig);
    //vec.push_back(&kEnergyScaleMuonInvSqrtUBBig);
    vec.push_back(&kEnergyScaleMuonFDBig);
    vec.push_back(&kEnergyScaleMuonSqrtFDBig);
    vec.push_back(&kEnergyScaleMuonInvSqrtFDBig);
    vec.push_back(&kEnergyScaleHadronBig);
    vec.push_back(&kEnergyScaleHadronSqrtBig);
    vec.push_back(&kEnergyScaleHadronInvSqrtBig);
    vec.push_back(&kEnergyScaleHadronNDBig);
    vec.push_back(&kEnergyScaleHadronSqrtNDBig);
    vec.push_back(&kEnergyScaleHadronInvSqrtNDBig);
    //vec.push_back(&kEnergyScaleHadronUBBig);
    //vec.push_back(&kEnergyScaleHadronSqrtUBBig);
    //vec.push_back(&kEnergyScaleHadronInvSqrtUBBig);
    vec.push_back(&kEnergyScaleHadronFDBig);
    vec.push_back(&kEnergyScaleHadronSqrtFDBig);
    vec.push_back(&kEnergyScaleHadronInvSqrtFDBig);

    return vec;
  }
} // namespace ana
