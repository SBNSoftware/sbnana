#include "sbnana/CAFAna/Core/FileReducer.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"

using namespace ana;

// Note the input file can (and likely should) be a SAM definition, or wildcard
// (if properly escaped from the shell).
void reduce_icarus(const std::string& in, const std::string& out)
{
  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSel = kNuMuCC_FullSelection;

  FileReducer reducer(in, out);
  reducer.AddReductionStep(ClearTrueParticles);
  reducer.AddSpillCut(kNumuSpillSel);
  reducer.AddSliceCut(kNumuSel);

  reducer.Go();
}
