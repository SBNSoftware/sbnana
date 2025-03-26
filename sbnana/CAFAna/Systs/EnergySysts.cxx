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

  void EnergyScaleSyst::Shift(double sigma, caf::SRTrueInteractionProxy *nu, double& weight) const
  {
  }

  void RecoEnergyScaleSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    //auto& det = sr->truth.det;
    bool all = (detector == EnergyScaleSystDetector::kAll);
    // We'd like to use these, but sr->truth.det doesn't seem to be working correctly right now
    //bool nd = (detector == EnergyScaleSystDetector::kSBND && det == caf::kSBND);
    //bool ub = (detector == EnergyScaleSystDetector::kMicroBooNE && det != caf::kSBND && det != caf::kICARUS);
    //bool fd = (detector == EnergyScaleSystDetector::kICARUS && det == caf::kICARUS);
    bool nd = (detector == EnergyScaleSystDetector::kSBND && sr->truth.baseline < 120);
    bool ub = (detector == EnergyScaleSystDetector::kMicroBooNE && sr->truth.baseline > 120 && sr->truth.baseline < 500);
    bool fd = (detector == EnergyScaleSystDetector::kSBND && sr->truth.baseline > 500);
    bool detector_cut = all || nd || ub || fd;
    if(!sr->truth.iscc || abs(sr->truth.pdg) != 14 || isnan(sr->fake_reco.nuE) || !detector_cut) 
      return ;

    for(auto& pfp: sr->reco.pfp) {
      if(pfp.trackScore < 0.5) continue;
      auto& trk = pfp.trk;
      if(part == EnergyScaleSystParticle::kMuon) {
        double scaleRange = sigma*uncertainty;
        double scaleMCS = sigma*uncertainty;
        switch(term) {
        case EnergyScaleSystTerm::kConstant:
          break;
        case EnergyScaleSystTerm::kSqrt:
          scaleRange *= std::sqrt(trk.rangeP.p_muon);
          if(!std::isnan(trk.mcsP.fwdP_muon))
            scaleMCS *= std::sqrt(trk.mcsP.fwdP_muon);
          break;
        case EnergyScaleSystTerm::kInverseSqrt:
          scaleRange /= std::sqrt(trk.rangeP.p_muon + 0.1);
          if(!std::isnan(trk.mcsP.fwdP_muon))
            scaleMCS /= std::sqrt(trk.mcsP.fwdP_muon + 0.1);
          break;
        }
        trk.rangeP.p_muon *= (1+scaleRange);
        if(!std::isnan(trk.mcsP.fwdP_muon))
          trk.mcsP.fwdP_muon *= (1+scaleMCS);
       } else if(part == EnergyScaleSystParticle::kHadron) {
        double scaleRangeP = sigma*uncertainty;
        double scaleMCSP = sigma*uncertainty;
        double scaleRangePi = sigma*uncertainty;
        double scaleMCSPi = sigma*uncertainty;
        switch(term) {
        case EnergyScaleSystTerm::kConstant:
          break;
        case EnergyScaleSystTerm::kSqrt:
          scaleRangeP *= std::sqrt(trk.rangeP.p_proton);
          if(!std::isnan(trk.mcsP.fwdP_proton))
            scaleMCSP *= std::sqrt(trk.mcsP.fwdP_proton);
          scaleRangePi *= std::sqrt(trk.rangeP.p_pion);
          if(!std::isnan(trk.mcsP.fwdP_pion))
            scaleMCSPi *= std::sqrt(trk.mcsP.fwdP_pion);
          break;
        case EnergyScaleSystTerm::kInverseSqrt:
          scaleRangeP /= std::sqrt(trk.rangeP.p_proton + 0.1);
          if(!std::isnan(trk.mcsP.fwdP_proton))
            scaleMCSP /= std::sqrt(trk.mcsP.fwdP_proton + 0.1);
          scaleRangePi /= std::sqrt(trk.rangeP.p_pion + 0.1);
          if(!std::isnan(trk.mcsP.fwdP_pion))
            scaleMCSPi /= std::sqrt(trk.mcsP.fwdP_pion + 0.1);
          break;
        }
        trk.rangeP.p_proton *= (1+scaleRangeP);
        if(!std::isnan(trk.mcsP.fwdP_proton))
          trk.mcsP.fwdP_proton *= (1+scaleMCSP);
        trk.rangeP.p_pion *= (1+scaleRangePi);
        if(!std::isnan(trk.mcsP.fwdP_pion))
          trk.mcsP.fwdP_pion *= (1+scaleMCSPi);
      }
    }
  }
  void RecoEnergyScaleSyst::Shift(double sigma, caf::SRTrueInteractionProxy *nu, double& weight) const
  {
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

   const RecoEnergyScaleSyst kRecoEnergyScaleMuon(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "RecoEnergyScaleMuon", "Correlated linear E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "RecoEnergyScaleMuonSqrt", "Correlated sqrt E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kAll, 0.02, "RecoEnergyScaleMuonInvSqrt", "Correlated inv sqrt E_{#mu} scale");

   const RecoEnergyScaleSyst kRecoEnergyScaleMuonND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "RecoEnergyScaleMuonND", "Uncorrelated SBND linear E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "RecoEnergyScaleMuonSqrtND", "Uncorrelated SBND sqrt E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kSBND, 0.02, "RecoEnergyScaleMuonInvSqrtND", "Uncorrelated SBND inv sqrt E_{#mu} scale");

   const RecoEnergyScaleSyst kRecoEnergyScaleMuonUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "RecoEnergyScaleMuonUB", "Uncorrelated MicroBooNE linear E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "RecoEnergyScaleMuonSqrtUB", "Uncorrelated MicroBooNE sqrt E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kMicroBooNE, 0.02, "RecoEnergyScaleMuonInvSqrtUB", "Uncorrelated MicroBooNE inv sqrt E_{#mu} scale");

   const RecoEnergyScaleSyst kRecoEnergyScaleMuonFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "RecoEnergyScaleMuonFD", "Uncorrelated ICARUS linear E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "RecoEnergyScaleMuonSqrtFD", "Uncorrelated ICARUS sqrt E_{#mu} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleMuonInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kMuon, EnergyScaleSystDetector::kICARUS, 0.02, "RecoEnergyScaleMuonInvSqrtFD", "Uncorrelated ICARUS inv sqrt E_{#mu} scale");

   const RecoEnergyScaleSyst kRecoEnergyScaleHadron(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "RecoEnergyScaleHadron", "Correlated linear E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrt(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "RecoEnergyScaleHadronSqrt", "Correlated sqrt E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrt(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kAll, 0.05, "RecoEnergyScaleHadronInvSqrt", "Correlated inv sqrt E_{had} scale");

   const RecoEnergyScaleSyst kRecoEnergyScaleHadronND(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "RecoEnergyScaleHadronND", "Uncorrelated SBND linear E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrtND(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "RecoEnergyScaleHadronSqrtND", "Uncorrelated SBND sqrt E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrtND(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kSBND, 0.05, "RecoEnergyScaleHadronInvSqrtND", "Uncorrelated SBND inv sqrt E_{had} scale");

   const RecoEnergyScaleSyst kRecoEnergyScaleHadronUB(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "RecoEnergyScaleHadronUB", "Uncorrelated MicroBooNE linear E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrtUB(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "RecoEnergyScaleHadronSqrtUB", "Uncorrelated MicroBooNE sqrt E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrtUB(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kMicroBooNE, 0.05, "RecoEnergyScaleHadronInvSqrtUB", "Uncorrelated MicroBooNE inv sqrt E_{had} scale");

   const RecoEnergyScaleSyst kRecoEnergyScaleHadronFD(EnergyScaleSystTerm::kConstant, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "RecoEnergyScaleHadronFD", "Uncorrelated ICARUS linear E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronSqrtFD(EnergyScaleSystTerm::kSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "RecoEnergyScaleHadronSqrtFD", "Uncorrelated ICARUS sqrt E_{had} scale");
   const RecoEnergyScaleSyst kRecoEnergyScaleHadronInvSqrtFD(EnergyScaleSystTerm::kInverseSqrt, EnergyScaleSystParticle::kHadron, EnergyScaleSystDetector::kICARUS, 0.05, "RecoEnergyScaleHadronInvSqrtFD", "Uncorrelated ICARUS inv sqrt E_{had} scale");

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

  std::vector<const ISyst*> GetRecoEnergySysts() {
    return {&kRecoEnergyScaleMuon,
            &kRecoEnergyScaleMuonSqrt,
            &kRecoEnergyScaleMuonInvSqrt,
            &kRecoEnergyScaleMuonND,
            &kRecoEnergyScaleMuonSqrtND,
            &kRecoEnergyScaleMuonInvSqrtND,
            &kRecoEnergyScaleMuonFD,
            &kRecoEnergyScaleMuonSqrtFD,
            &kRecoEnergyScaleMuonInvSqrtFD,
            &kRecoEnergyScaleHadron,
            &kRecoEnergyScaleHadronSqrt,
            &kRecoEnergyScaleHadronInvSqrt,
            &kRecoEnergyScaleHadronND,
            &kRecoEnergyScaleHadronSqrtND,
            &kRecoEnergyScaleHadronInvSqrtND,
            &kRecoEnergyScaleHadronFD,
            &kRecoEnergyScaleHadronSqrtFD,
            &kRecoEnergyScaleHadronInvSqrtFD};
  }
} // namespace ana
