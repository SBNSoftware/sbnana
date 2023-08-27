#pragma once

#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecRecoEventVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecTrueEventVars.h"

using namespace ana;
using namespace std;

namespace ana{

  std::vector<std::string> GetNuMITrueTreeLabels();
  std::vector<SpillMultiVar> GetNuMITrueTreeVars();

  std::vector<std::string> GetNuMIRecoTreeLabels();
  std::vector<Var> GetNuMIRecoTreeVars();

}
