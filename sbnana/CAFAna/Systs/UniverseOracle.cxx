#include "sbnana/CAFAna/Systs/UniverseOracle.h"

#include "sbnanaobj/StandardRecord/SRGlobal.h"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSeqCollection.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  // Go grubbing through all the open files looking for a globalTree to pull an
  // SRGlobal object from. This is a hack, probably SpectrumLoader should have
  // some hooks to send SRGlobal to interested parties each time it opens a new
  // file(?)
  caf::SRGlobal GetSRGlobal()
  {
    caf::SRGlobal global;
    bool got = false;

    TSeqCollection* seq = gROOT->GetListOfFiles();
    for(int i = 0; i < seq->GetEntries(); ++i){
      TFile* f = (TFile*)seq->At(i);
      if(f->GetOption() != std::string("READ")) continue;
      TTree* tr = (TTree*)f->Get("globalTree");
      if(!tr) continue;
      if(tr->GetEntries() < 1) continue;

      if(got){
        std::cout << "\nUniverseOracle: Found globalTree in multiple places. Will use first one, but this is unexpected" << std::endl;
      }

      caf::SRGlobal* pglobal = &global;
      tr->SetBranchAddress("global", &pglobal);
      tr->GetEntry(0);
      got = true;
    }

    if(!got){
      std::cout << "\nUniverseOracle: Failed to find globalTree in any open input file" << std::endl;
      abort();
    }

    return global;
  }

  // --------------------------------------------------------------------------
  void PrintSRGlobal(const caf::SRGlobal& global)
  {
    std::cout << global.wgts.size() << " parameter sets:" << std::endl;
    for(unsigned int i = 0; i < global.wgts.size(); ++i){
      const caf::SRWeightPSet& pset = global.wgts[i];
      std::cout << "  " << i << ": " << pset.name << ", type " << pset.type << ", " << pset.nuniv << " universes, adjusted parameters:" << std::endl;

      for(const caf::SRWeightMapEntry& entry: pset.map){
        std::cout << "    " << entry.param.name << std::endl;
      }
    }
  }

  // --------------------------------------------------------------------------
  UniverseOracle& UniverseOracle::Instance()
  {
    static UniverseOracle uo;
    return uo;
  }

  // --------------------------------------------------------------------------
  UniverseOracle::UniverseOracle()
  {
    const caf::SRGlobal global = GetSRGlobal();
    std::cout << "\nSystematic weights in file:" << std::endl;
    PrintSRGlobal(global);

    for(unsigned int i = 0; i < global.wgts.size(); ++i){
      const caf::SRWeightPSet& pset = global.wgts[i];

      // Save the pset index in all cases
      fPSetIdxs[pset.name] = i;

      // Only save the remaining fields in parameter sets where only a single
      // knob is shifted
      if(pset.map.size() != 1) continue;

      // Save which position in the vector this was
      fSystIdxs[pset.name] = i;
      // Save all the knob values
      fShiftVals[pset.name] = pset.map[0].vals;
    }
  }

  // --------------------------------------------------------------------------
  bool UniverseOracle::SystExists(const std::string& name) const
  {
    return fShiftVals.find(name) != fShiftVals.end();
  }

  // --------------------------------------------------------------------------
  std::vector<std::string> UniverseOracle::Systs() const
  {
    std::vector<std::string> ret;
    ret.reserve(fShiftVals.size());
    for(auto it: fShiftVals) ret.push_back(it.first);
    return ret;
  }

  // --------------------------------------------------------------------------
  const std::vector<float>& UniverseOracle::ShiftsForSyst(const std::string& name) const
  {
    assert(SystExists(name));
    return fShiftVals.find(name)->second;
  }

  // --------------------------------------------------------------------------
  unsigned int UniverseOracle::ParameterSetIndex(const std::string& name) const
  {
    auto it = fPSetIdxs.find(name);
    if(it == fPSetIdxs.end()){
      std::cout << "UniverseOracle: pset '" << name << "' not known" << std::endl;
      abort();
    }
    return it->second;
  }

  // --------------------------------------------------------------------------
  unsigned int UniverseOracle::SystIndex(const std::string& name) const
  {
    auto it = fSystIdxs.find(name);
    if(it == fSystIdxs.end()){
      std::cout << "UniverseOracle: syst '" << name << "' not known" << std::endl;
      abort();
    }
    return it->second;
  }

  // --------------------------------------------------------------------------
  unsigned int UniverseOracle::ClosestShiftIndex(const std::string& name,
                                                 double shift,
                                                 ESide side,
                                                 double* trueShift) const
  {
    const std::vector<float>& v = ShiftsForSyst(name);
    int bestIdx = -1;
    double bestDist;
    for(unsigned int i = 0; i < v.size(); ++i){
      if(side == ESide::kBelow && v[i] > shift) continue;
      if(side == ESide::kAbove && v[i] < shift) continue;
      const double dv = fabs(v[i]-shift);
      if(bestIdx == -1 || dv < bestDist){
        bestIdx = i;
        bestDist = dv;
        if(trueShift) *trueShift = v[i];
      }
    }
    return unsigned(bestIdx);
  }
}
