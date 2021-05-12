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

  // Muon energy scale syst (correlated among all three experiments)
  class EnergyScaleMuon: public ISyst
  {
  public:
  EnergyScaleMuon() : ISyst("EnergyScaleMuon", "2% Scale Uncertainty on Muon Energy") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.02 * sigma;
      //For dummying up from existing info, use truth
      //if (numu cc) {
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && !isnan(sr->freco.lepton.ke)) {
	//reco_energy = (reco_energy - muon_energy) + muon_energy(1. + scale);
	sr->freco.nuE += sr->freco.lepton.ke * scale;
      }
    }
  };
  extern const EnergyScaleMuon kEnergyScaleMuon;

  //////////////////////////////////////////////

  // Muon energy scale syst (SBND only)
  class EnergyScaleMuonND: public ISyst
  {
  public:
  EnergyScaleMuonND() : ISyst("EnergyScaleMuonND", "2% Scale Uncertainty on Muon Energy in SBND") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.02 * sigma;
      //For dummying up from existing info, use truth
      //if (numu cc) {
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && !isnan(sr->freco.lepton.ke) && int(sr->truth.baseline) < 120) {
	//reco_energy = (reco_energy - muon_energy) + muon_energy(1. + scale);
	sr->freco.nuE += sr->freco.lepton.ke * scale;
      }
    }
  };
  extern const EnergyScaleMuonND kEnergyScaleMuonND;

  //////////////////////////////////////////////

//  // Muon energy scale syst (uBooNE only)
//  class EnergyScaleMuonMB: public ISyst
//  {
//  public:
//  EnergyScaleMuonMB() : ISyst("EnergyScaleMuonMB", "2% Scale Uncertainty on Muon Energy in uBooNE") {}
//    void Shift(double sigma,
//	       caf::SRSliceProxy* sr, double& weight) const override
//    {
//      double scale = 0.02 * sigma;
//      //For dummying up from existing info, use truth
//      //if (numu cc) {
//      if (sr->truth[0].neutrino.iscc and abs(sr->truth[0].neutrino.pdg)==14 and int(sr->experiment) == kMicroBoone) {
//	//reco_energy = (reco_energy - muon_energy) + muon_energy(1. + scale);
//	sr->reco.reco_energy += sr->truth[0].lepton.energy * scale;
//      }
//    }
//  };
//  extern const EnergyScaleMuonMB kEnergyScaleMuonMB;
//
  //////////////////////////////////////////////

  // Muon energy scale syst (ICARUS only)
  class EnergyScaleMuonFD: public ISyst
  {
  public:
  EnergyScaleMuonFD() : ISyst("EnergyScaleMuonFD", "2% Scale Uncertainty on Muon Energy in ICARUS") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.02 * sigma;
      //For dummying up from existing info, use truth
      //if (numu cc) {
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && int(sr->truth.baseline) > 500 && !isnan(sr->freco.lepton.ke)) {
	//reco_energy = (reco_energy - muon_energy) + muon_energy(1. + scale);
	sr->freco.nuE += sr->freco.lepton.ke * scale;
      }
    }
  };
  extern const EnergyScaleMuonFD kEnergyScaleMuonFD;

  // Muon energy scale syst (ICARUS only)
  class EnergyScaleMuonFD5p: public ISyst
  {
  public:
  EnergyScaleMuonFD5p() : ISyst("EnergyScaleMuonFD5p", "5% Scale Uncertainty on Muon Energy in ICARUS") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.05 * sigma;
      //For dummying up from existing info, use truth
      //if (numu cc) {
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && int(sr->truth.baseline) > 500 && !isnan(sr->freco.lepton.ke)) {
	//reco_energy = (reco_energy - muon_energy) + muon_energy(1. + scale);
	sr->freco.nuE += sr->freco.lepton.ke * scale;
      }
    }
  };
  extern const EnergyScaleMuonFD5p kEnergyScaleMuonFD5p;

  // Muon energy scale syst (ICARUS only)
  class EnergyScaleMuonFD10p: public ISyst
  {
  public:
  EnergyScaleMuonFD10p() : ISyst("EnergyScaleMuonFD10p", "10% Scale Uncertainty on Muon Energy in ICARUS") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.10 * sigma;
      //For dummying up from existing info, use truth
      //if (numu cc) {
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && int(sr->truth.baseline) > 500 && !isnan(sr->freco.lepton.ke)) {
	//reco_energy = (reco_energy - muon_energy) + muon_energy(1. + scale);
	sr->freco.nuE += sr->freco.lepton.ke * scale;
      }
    }
  };
  extern const EnergyScaleMuonFD10p kEnergyScaleMuonFD10p;

  // Muon energy scale syst (ICARUS only)
  class EnergyScaleMuonFD20p: public ISyst
  {
  public:
  EnergyScaleMuonFD20p() : ISyst("EnergyScaleMuonFD20p", "20% Scale Uncertainty on Muon Energy in ICARUS") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.20 * sigma;
      //For dummying up from existing info, use truth
      //if (numu cc) {
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && int(sr->truth.baseline) > 500 && !isnan(sr->freco.lepton.ke)) {
	//reco_energy = (reco_energy - muon_energy) + muon_energy(1. + scale);
	sr->freco.nuE += sr->freco.lepton.ke * scale;
      }
    }
  };
  extern const EnergyScaleMuonFD20p kEnergyScaleMuonFD20p;
  //////////////////////////////////////////////

  // Hadron energy scale syst (correlated among all three experiments)
  class EnergyScaleHadron: public ISyst
  {
  public:
  EnergyScaleHadron() : ISyst("EnergyScaleHadron", "5% Scale Uncertainty on Hadron Energy") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.05 * sigma;
      //For dummying up from existing info, use truth
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && sr->freco.nhad > 0) {
        for (size_t i = 0; i < sr->freco.hadrons.size(); ++i) {
          sr->freco.nuE += sr->freco.hadrons[i].ke * scale;
        }
	//sr->reco.reco_energy += (sr->reco.reco_energy - sr->truth[0].lepton.energy) * scale;
      }
    }
  };
  extern const EnergyScaleHadron kEnergyScaleHadron;

  //////////////////////////////////////////////

  // Hadron energy scale syst (SBND only)
  class EnergyScaleHadronND: public ISyst
  {
  public:
  EnergyScaleHadronND() : ISyst("EnergyScaleHadronND", "5% Scale Uncertainty on Hadron Energy in SBND") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.05 * sigma;
      //For dummying up from existing info, use truth
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && sr->freco.nhad > 0 && int(sr->truth.baseline) < 120) {
        for (size_t i = 0; i < sr->freco.hadrons.size(); ++i) {
          sr->freco.nuE += sr->freco.hadrons[i].ke * scale;
        }
	//sr->reco.reco_energy += (sr->reco.reco_energy - sr->truth[0].lepton.energy) * scale;
      }
    }
  };
  extern const EnergyScaleHadronND kEnergyScaleHadronND;

//  //////////////////////////////////////////////
//
//  // Hadron energy scale syst (uBOONE only)
//  class EnergyScaleHadronMB: public ISyst
//  {
//  public:
//  EnergyScaleHadronMB() : ISyst("EnergyScaleHadronMB", "5% Scale Uncertainty on Hadron Energy in uBOONE") {}
//    void Shift(double sigma,
//	       caf::SRSliceProxy* sr, double& weight) const override
//    {
//      double scale = 0.05 * sigma;
//      //For dummying up from existing info, use truth
//      if (sr->truth[0].neutrino.iscc and abs(sr->truth[0].neutrino.pdg)==14 and int(sr->experiment) == kMicroBoone) {
//	sr->reco.reco_energy += (sr->reco.reco_energy - sr->truth[0].lepton.energy) * scale;
//      }
//    }
//  };
//  extern const EnergyScaleHadronMB kEnergyScaleHadronMB;
//
//  //////////////////////////////////////////////
//
  // Hadron energy scale syst (ICARUS only)
  class EnergyScaleHadronFD: public ISyst
  {
  public:
  EnergyScaleHadronFD() : ISyst("EnergyScaleHadronFD", "5% Scale Uncertainty on Hadron Energy in ICARUS") {}
    void Shift(double sigma,
	       caf::SRSliceProxy* sr, double& weight) const override
    {
      double scale = 0.05 * sigma;
      //For dummying up from existing info, use truth
      if (sr->truth.iscc && abs(sr->truth.pdg)==14 && sr->freco.nhad > 0 && int(sr->truth.baseline) > 500) {
        for (size_t i = 0; i < sr->freco.hadrons.size(); ++i) {
          sr->freco.nuE += sr->freco.hadrons[i].ke * scale;
        }
	//sr->reco.reco_energy += (sr->reco.reco_energy - sr->truth[0].lepton.energy) * scale;
      }
    }
  };
  extern const EnergyScaleHadronFD kEnergyScaleHadronFD;

  //////////////////////////////////////////////

  // Vector of energy scale systematics
  struct EnergySystVector: public std::vector<const ISyst*>
  {

  };

EnergySystVector GetEnergySysts();

}
