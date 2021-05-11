#include "sbnana/CAFAna/Core/FileReducer.h"

#include "sbnana/SBNAna/Cuts/NumuCutsSBND202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

using namespace ana;

// Since I'm in the SBND group the SAM dataset works fine for me
//
// NB this has to be the NON FLAT version for FileReducer to work
const std::string sbnd_wildcard = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_nu_spill_gsimple-configf-v1_tpc_caf_sbnd";

void reduce_sbnd()
{
  const SpillCut kNumuSpillSel = kNoSpillCut; // TODO
  const Cut kNumuSel = kSlcIsRecoNu && kSlcNuScoreCut && kFiducialVolumeND && kSlcFlashMatchCut && kHasPrimaryMuonTrk && kCRTTrackAngleCut && kCRTHitDistanceCut;

  FileReducer reducer(sbnd_wildcard, "reduced_sbnd.root");
  reducer.AddSpillCut(kNumuSpillSel);
  reducer.AddSliceCut(kNumuSel);

  reducer.Go();
}
