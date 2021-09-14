#include "sbnana/CAFAna/Systs/EnergySysts.h"

namespace ana {
   const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuon(0.02, "EnergyScaleMuon");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonSqrt(0.02, "EnergyScaleMuonSqrt");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonInvSqrt(0.02, "EnergyScaleMuonInvSqrt");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonND(0.02, "EnergyScaleMuonND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtND(0.02, "EnergyScaleMuonSqrtND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtND(0.02, "EnergyScaleMuonInvSqrtND");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonFD(0.02, "EnergyScaleMuonFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtFD(0.02, "EnergyScaleMuonSqrtFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtFD(0.02, "EnergyScaleMuonInvSqrtFD");

   const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadron(0.05, "EnergyScaleHadron");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronSqrt(0.05, "EnergyScaleHadronSqrt");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronInvSqrt(0.05, "EnergyScaleHadronInvSqrt");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronND(0.05, "EnergyScaleHadronND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtND(0.05, "EnergyScaleHadronSqrtND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtND(0.05, "EnergyScaleHadronInvSqrtND");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronFD(0.05, "EnergyScaleHadronFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtFD(0.05, "EnergyScaleHadronSqrtFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtFD(0.05, "EnergyScaleHadronInvSqrtFD");

   const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonBig(0.05, "EnergyScaleMuonBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonSqrtBig(0.05, "EnergyScaleMuonSqrtBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Muon>        kEnergyScaleMuonInvSqrtBig(0.05, "EnergyScaleMuonInvSqrtBig");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonNDBig(0.05, "EnergyScaleMuonNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtNDBig(0.05, "EnergyScaleMuonSqrtNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtNDBig(0.05, "EnergyScaleMuonInvSqrtNDBig");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonFDBig(0.05, "EnergyScaleMuonFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonSqrtFDBig(0.05, "EnergyScaleMuonSqrtFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Muon>    kEnergyScaleMuonInvSqrtFDBig(0.05, "EnergyScaleMuonInvSqrtFDBig");

   const EnergyScaleSystCorr<EnergyScaleSystTerm::Constant, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronBig(0.10, "EnergyScaleHadronBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::Sqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronSqrtBig(0.10, "EnergyScaleHadronSqrtBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::InverseSqrt, 
                                   EnergyScaleSystParticle::Hadron>      kEnergyScaleHadronInvSqrtBig(0.10, "EnergyScaleHadronInvSqrtBig");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronNDBig(0.10, "EnergyScaleHadronNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtNDBig(0.10, "EnergyScaleHadronSqrtNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtNDBig(0.10, "EnergyScaleHadronInvSqrtNDBig");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Constant, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronFDBig(0.10, "EnergyScaleHadronFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::Sqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronSqrtFDBig(0.10, "EnergyScaleHadronSqrtFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::InverseSqrt, 
                                       EnergyScaleSystParticle::Hadron>  kEnergyScaleHadronInvSqrtFDBig(0.10, "EnergyScaleHadronInvSqrtFDBig");

  EnergySystVector GetEnergySysts() {

    EnergySystVector vec;
    vec.push_back(&kEnergyScaleMuon);
    vec.push_back(&kEnergyScaleMuonSqrt);
    vec.push_back(&kEnergyScaleMuonInvSqrt);
    vec.push_back(&kEnergyScaleMuonND);
    vec.push_back(&kEnergyScaleMuonSqrtND);
    vec.push_back(&kEnergyScaleMuonInvSqrtND);
    vec.push_back(&kEnergyScaleMuonFD);
    vec.push_back(&kEnergyScaleMuonSqrtFD);
    vec.push_back(&kEnergyScaleMuonInvSqrtFD);
    vec.push_back(&kEnergyScaleHadron);
    //vec.push_back(&kEnergyScaleChargedHadron);
    //vec.push_back(&kEnergyScalePi0);
    //vec.push_back(&kEnergyScaleNeutron);
    vec.push_back(&kEnergyScaleHadronSqrt);
    //vec.push_back(&kEnergyScalePi0Sqrt);
    //vec.push_back(&kEnergyScaleNeutronSqrt);
    vec.push_back(&kEnergyScaleHadronInvSqrt);
    //vec.push_back(&kEnergyScalePi0InvSqrt);
    //vec.push_back(&kEnergyScaleNeutronInvSqrt);
    vec.push_back(&kEnergyScaleHadronND);
    //vec.push_back(&kEnergyScaleChargedHadronND);
    //vec.push_back(&kEnergyScalePi0ND);
    //vec.push_back(&kEnergyScaleNeutronND);
    vec.push_back(&kEnergyScaleHadronSqrtND);
    //vec.push_back(&kEnergyScalePi0SqrtND);
    //vec.push_back(&kEnergyScaleNeutronSqrtND);
    vec.push_back(&kEnergyScaleHadronInvSqrtND);
    //vec.push_back(&kEnergyScalePi0InvSqrtND);
    //vec.push_back(&kEnergyScaleNeutronInvSqrtND);
    vec.push_back(&kEnergyScaleHadronFD);
    //vec.push_back(&kEnergyScaleChargedHadronFD);
    //vec.push_back(&kEnergyScalePi0FD);
    //vec.push_back(&kEnergyScaleNeutronFD);
    vec.push_back(&kEnergyScaleHadronSqrtFD);
    //vec.push_back(&kEnergyScalePi0SqrtFD);
    //vec.push_back(&kEnergyScaleNeutronSqrtFD);
    vec.push_back(&kEnergyScaleHadronInvSqrtFD);
    //vec.push_back(&kEnergyScalePi0InvSqrtFD);
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
    vec.push_back(&kEnergyScaleMuonFDBig);
    vec.push_back(&kEnergyScaleMuonSqrtFDBig);
    vec.push_back(&kEnergyScaleMuonInvSqrtFDBig);
    vec.push_back(&kEnergyScaleHadronBig);
    vec.push_back(&kEnergyScaleHadronSqrtBig);
    vec.push_back(&kEnergyScaleHadronInvSqrtBig);
    vec.push_back(&kEnergyScaleHadronNDBig);
    vec.push_back(&kEnergyScaleHadronSqrtNDBig);
    vec.push_back(&kEnergyScaleHadronInvSqrtNDBig);
    vec.push_back(&kEnergyScaleHadronFDBig);
    vec.push_back(&kEnergyScaleHadronSqrtFDBig);
    vec.push_back(&kEnergyScaleHadronInvSqrtFDBig);

    return vec;
  }
} // namespace ana
