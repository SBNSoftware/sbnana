#include "sbnana/CAFAna/Core/FileReducer.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"

using namespace ana;

// This should be
// workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus
// but that is only available from SAM to icarus users, for now.
const std::string icarus_wildcard = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/*/*/*.root";

void reduce_icarus()
{
  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSel = kNuMuCC_FullSelection;

  FileReducer reducer(icarus_wildcard, "reduced_icarus.root");
  reducer.AddSpillCut(kNumuSpillSel);
  reducer.AddSliceCut(kNumuSel);

  reducer.Go();
}
