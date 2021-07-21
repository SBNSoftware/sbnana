#pragma once

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Ratio.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/NueCuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"
#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Vars/Vars.h"

#include <fstream>

using namespace ana;

const SpillVar kNuEnergy([](const caf::SRSpillProxy* sr) {
  float energy = (sr->mc.nnu != 1 ? -5.f : (float)sr->mc.nu[0].E);
  return energy;
});

const SpillCut kSpillSingleNu([](const caf::SRSpillProxy* sr) {
  return (sr->mc.nnu < 2);
});

const Cut kTrueNuFV(
    [](const caf::SRSliceProxy* slc) {
      if (slc->truth.index == -1)
        return false;

      return PtInVolAbsX(slc->truth.position, fvndAbs);
    });

const Cut kNuEEnergyCut(
    [](const caf::SRSliceProxy* slc) {
      if (slc->truth.index == -1)
        return true;

      for (auto const& prim : slc->truth.prim) {
        if (std::abs(prim.pdg) != 11 && std::abs(prim.pdg) != 22)
          continue;

        if (prim.startE > 0.2)
          return true;
      }

      return false;
    });

const Cut kCosmicRay = !kHasNu;

// In this example, our signal is Nue cc
const Cut kNuECC = kIsNue && !kIsNC && kSlcCompletenessCut && kTrueNuFV && kNuEEnergyCut;
const Cut kNuMuCC = kIsNumu && !kIsNC && kSlcCompletenessCut && kTrueNuFV;
const Cut kNC = kIsNC && kSlcCompletenessCut && kTrueNuFV;
const Cut kOtherNu = kHasNu && !(kNuECC || kNuMuCC || kNC);

double POT = 6.6E20;
