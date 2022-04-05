#include "sbnana/CAFAna/Systs/EnergySysts.h"

namespace ana {
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuon(0.02, "EnergyScaleMuon");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonSqrt(0.02, "EnergyScaleMuonSqrt");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonInvSqrt(0.02, "EnergyScaleMuonInvSqrt");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonND(0.02, "EnergyScaleMuonND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtND(0.02, "EnergyScaleMuonSqrtND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtND(0.02, "EnergyScaleMuonInvSqrtND");

   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonUB(0.02, "EnergyScaleMuonUB");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtUB(0.02, "EnergyScaleMuonSqrtUB");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtUB(0.02, "EnergyScaleMuonInvSqrtUB");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonFD(0.02, "EnergyScaleMuonFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtFD(0.02, "EnergyScaleMuonSqrtFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtFD(0.02, "EnergyScaleMuonInvSqrtFD");

   const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadron(0.05, "EnergyScaleHadron");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronSqrt(0.05, "EnergyScaleHadronSqrt");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronInvSqrt(0.05, "EnergyScaleHadronInvSqrt");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronND(0.05, "EnergyScaleHadronND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtND(0.05, "EnergyScaleHadronSqrtND");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtND(0.05, "EnergyScaleHadronInvSqrtND");

   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronUB(0.05, "EnergyScaleHadronUB");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtUB(0.05, "EnergyScaleHadronSqrtUB");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtUB(0.05, "EnergyScaleHadronInvSqrtUB");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronFD(0.05, "EnergyScaleHadronFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtFD(0.05, "EnergyScaleHadronSqrtFD");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtFD(0.05, "EnergyScaleHadronInvSqrtFD");

   const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonBig(0.05, "EnergyScaleMuonBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonSqrtBig(0.05, "EnergyScaleMuonSqrtBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kMuon>        kEnergyScaleMuonInvSqrtBig(0.05, "EnergyScaleMuonInvSqrtBig");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonNDBig(0.05, "EnergyScaleMuonNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtNDBig(0.05, "EnergyScaleMuonSqrtNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtNDBig(0.05, "EnergyScaleMuonInvSqrtNDBig");

   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonUBBig(0.05, "EnergyScaleMuonUBBig");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtUBBig(0.05, "EnergyScaleMuonSqrtUBBig");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtUBBig(0.05, "EnergyScaleMuonInvSqrtUBBig");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonFDBig(0.05, "EnergyScaleMuonFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonSqrtFDBig(0.05, "EnergyScaleMuonSqrtFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kMuon>    kEnergyScaleMuonInvSqrtFDBig(0.05, "EnergyScaleMuonInvSqrtFDBig");

   const EnergyScaleSystCorr<EnergyScaleSystTerm::kConstant, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronBig(0.10, "EnergyScaleHadronBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronSqrtBig(0.10, "EnergyScaleHadronSqrtBig");
   const EnergyScaleSystCorr<EnergyScaleSystTerm::kInverseSqrt, 
                                   EnergyScaleSystParticle::kHadron>      kEnergyScaleHadronInvSqrtBig(0.10, "EnergyScaleHadronInvSqrtBig");

   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronNDBig(0.10, "EnergyScaleHadronNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtNDBig(0.10, "EnergyScaleHadronSqrtNDBig");
   const EnergyScaleSystUncorrND<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtNDBig(0.10, "EnergyScaleHadronInvSqrtNDBig");

   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronUBBig(0.10, "EnergyScaleHadronUBBig");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtUBBig(0.10, "EnergyScaleHadronSqrtUBBig");
   const EnergyScaleSystUncorrUB<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtUBBig(0.10, "EnergyScaleHadronInvSqrtUBBig");

   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kConstant, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronFDBig(0.10, "EnergyScaleHadronFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronSqrtFDBig(0.10, "EnergyScaleHadronSqrtFDBig");
   const EnergyScaleSystUncorrFD<EnergyScaleSystTerm::kInverseSqrt, 
                                       EnergyScaleSystParticle::kHadron>  kEnergyScaleHadronInvSqrtFDBig(0.10, "EnergyScaleHadronInvSqrtFDBig");

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
