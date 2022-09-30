#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202208.h"
#include "sbnana/SBNAna/Vars/NumuVarsIcarus202208.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <fstream>

using namespace ana;

void selected_events()
{
  SpectrumLoader loader("icarus_BNB_Nu_Cosmics_v09_37_02_04_flatcaf");

  std::ofstream out("selected_events.txt");

  const SpillVar dummy_var([&out](const caf::SRSpillProxy* spill){
    for(const auto& slc: spill->slc) {
      if(kIcarus202208QELikeContainedMuon(&slc)) {
        out << "Run: " << spill->hdr.run << "\tSubrun: " << spill->hdr.subrun
            << "\tEvent: " << spill->hdr.evt << "\tVertex: (" 
            << slc.vertex.x << ", " << slc.vertex.y << ", " << slc.vertex.z << ")\n";
        return 1;
      }
    } 
    return 0;
  });

  const Binning dummy_bins = Binning::Simple(2, 0, 2);
  Spectrum dummy_spec("Dummy Label", dummy_bins, loader, dummy_var, kNoSpillCut);

  loader.Go();
}
