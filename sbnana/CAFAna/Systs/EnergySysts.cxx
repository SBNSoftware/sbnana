#include "sbnana/CAFAna/Systs/EnergySysts.h"

namespace ana {
  const EnergyScaleMuon kEnergyScaleMuon;
  const EnergyScaleMuonND kEnergyScaleMuonND;
  // const EnergyScaleMuonMB kEnergyScaleMuonMB;
   const EnergyScaleMuonFD kEnergyScaleMuonFD;
   const EnergyScaleMuonFD5p kEnergyScaleMuonFD5p;
   const EnergyScaleMuonFD10p kEnergyScaleMuonFD10p;
   const EnergyScaleMuonFD20p kEnergyScaleMuonFD20p;
   const EnergyScaleHadron kEnergyScaleHadron;
   const EnergyScaleHadronND kEnergyScaleHadronND;
  // const EnergyScaleHadronMB kEnergyScaleHadronMB;
   const EnergyScaleHadronFD kEnergyScaleHadronFD;
  // const UncorrFDTotSqrt kUncorrFDTotSqrt;
  // const UncorrFDTotInvSqrt kUncorrFDTotInvSqrt;
  // const UncorrFDHadSqrt kUncorrFDHadSqrt;
  // const UncorrFDHadInvSqrt kUncorrFDHadInvSqrt;
  // const UncorrFDMuSqrt kUncorrFDMuSqrt;
  // const UncorrFDMuInvSqrt kUncorrFDMuInvSqrt;
  // const UncorrFDNSqrt kUncorrFDNSqrt;
  // const UncorrFDNInvSqrt kUncorrFDNInvSqrt;
  // const UncorrFDEMSqrt kUncorrFDEMSqrt;
  // const UncorrFDEMInvSqrt kUncorrFDEMInvSqrt;
  // const EScaleMuLArFD kEScaleMuLArFD;
  // const ChargedHadUncorrFD kChargedHadUncorrFD;
  // const NUncorrFD kNUncorrFD;
  // const EMUncorrFD kEMUncorrFD;
  // const MuonResFD kMuonResFD;
  // const EMResFD kEMResFD;
  // const ChargedHadResFD kChargedHadResFD;
  // const NResFD kNResFD;

  EnergySystVector GetEnergySysts() {

    EnergySystVector vec;
    vec.push_back(&kEnergyScaleMuon);
    vec.push_back(&kEnergyScaleMuonND);
    // vec.push_back(&kEnergyScaleMuonMB);
    vec.push_back(&kEnergyScaleMuonFD);
    vec.push_back(&kEnergyScaleMuonFD5p);
    vec.push_back(&kEnergyScaleMuonFD10p);
    vec.push_back(&kEnergyScaleMuonFD20p);
    vec.push_back(&kEnergyScaleHadron);
    vec.push_back(&kEnergyScaleHadronND);
    // vec.push_back(&kEnergyScaleHadronMB);
    vec.push_back(&kEnergyScaleHadronFD);

    return vec;
  }
} // namespace ana
