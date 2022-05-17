const Cut kNonUnambiguousSlice([](const caf::SRSliceProxy *slc)
     {
       return !slc->is_clear_cosmic;
     });

const SpillCut kHasNonUnambiguousSlice([](const caf::SRSpillProxy *sr)
     {
       for(auto const &slc : sr->slc)
         {
           if(!slc.is_clear_cosmic) return true;
         }

       return false;
     });

const Var kCRUMBSScore = SIMPLEVAR(crumbs_result.score);

const SpillVar kBestCRUMBSScore([](const caf::SRSpillProxy *sr) -> double
     {
       double bestValue = -std::numeric_limits<double>::max();
       
       for(auto const &slc : sr->slc)
         {
           if(!slc.is_clear_cosmic)
             {
               if(slc.crumbs_result.score > bestValue)
                 bestValue = slc.crumbs_result.score;
             }
         }
       
       return bestValue;
     });

const SpillVar kBestCRUMBSSliceID([](const caf::SRSpillProxy *sr) -> double
     {
       unsigned id      = 0;
       unsigned bestId  = std::numeric_limits<unsigned>::max();
       double bestValue = -std::numeric_limits<double>::max();

       for(auto const &slc : sr->slc)
         {
           if(!slc.is_clear_cosmic)
             {
               if(slc.crumbs_result.score > bestValue)
                 {
                   bestValue = slc.crumbs_result.score;
                   bestId    = id;
                 }
             }
           ++id;
         }

       return bestId;
     });

const SpillVar kBestCRUMBSSliceNTracks([](const caf::SRSpillProxy *sr) -> float
     {
       return sr->slc[kBestCRUMBSSliceID(sr)].reco.ntrk;
     });

const Cut kCRUMBSCut([](const caf::SRSliceProxy *slc)
     {
       return slc->crumbs_result.score > -0.05;
     });

const SpillCut kBestCRUMBSSliceCut([](const caf::SRSpillProxy *sr)
     {
       return sr->slc[kBestCRUMBSSliceID(sr)].crumbs_result.score > -0.05;
     });
