#include "sbnana/SBNAna/Cuts/ICARUSDataQualityCuts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

namespace ana::icarus
{
  // Laser Runs
  const ana::SpillCut kLaserRunsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isLaserRun_run2(spill->hdr.run);
    }
  );

  // TPC Temp
  const ana::SpillCut kTPCTempRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not hasTPCTempIssue_run2(spill->hdr.run);
    }
  );

  // Trigger Tetss
  const ana::SpillCut kTriggerTestsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isTriggerTest_run2(spill->hdr.run);
    }
  );

  // Failed Runs
  const ana::SpillCut kFailedRunsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isFailedRun_run2(spill->hdr.run);
    }
  );

  // Unspecified Tests
  const ana::SpillCut kUnspecifiedTestsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isUnspecifiedTest_run2(spill->hdr.run);
    }
  );

  // Server Issues
  const ana::SpillCut kServerIssuesRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not hasServerIssues_run2(spill->hdr.run);
    }
  );

  // Tagged By RunCo
  const ana::SpillCut kTaggedByRunCoRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isTaggedByRunCo(spill->hdr.run);
    }
  );

  // PMT Tests
  const ana::SpillCut kPMTTestsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isPMTTest_run2(spill->hdr.run);
    }
  );

  // CRT Tests
  const ana::SpillCut kCRTTestsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isCRTTest_run2(spill->hdr.run);
    }
  );

  // Absent from RunCo Sheet
  const ana::SpillCut kAbsentFromRunCoSheetRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isNotInSheet_run2(spill->hdr.run);
    }
  );

  // DAQ Tests
  const ana::SpillCut kDAQTestsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isDAQTest_run2(spill->hdr.run);
    }
  );

  // TPC Tests
  const ana::SpillCut kTPCTestsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isTPCTest_run2(spill->hdr.run);
    }
  );

  // WireBias PS Problems
  const ana::SpillCut kWireBiasRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not hasWireBiasIssue_run2(spill->hdr.run);
    }
  );

  // Short Runs
  const ana::SpillCut kShortRunsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isShortRun_run2(spill->hdr.run);
    }
  );

  // Missing TPC Components
  const ana::SpillCut kMissingTPCRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isMissingTPC_run2(spill->hdr.run);
    }
  );

  // Invalid Timing Information
  const ana::SpillCut kInvalidTimingRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isMissingTimingInfo_run2(spill->hdr.run);
    }
  );

  // Absent from Database
  const ana::SpillCut kAbsentDBRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isNotInDatabase_run2(spill->hdr.run);
    }
  );

  // Test Configuration (Test_thr390_LockTemp_true_00001 only)
  const ana::SpillCut kTestConfigRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isTestConfig_run2(spill->hdr.run);
    }
  );

  // Missing PMT Components
  const ana::SpillCut kMissingPMTRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not isMissingPMT_run2(spill->hdr.run);
    }
  );

  // Prefilter
  const ana::SpillCut kPrefilterRun2 = kLaserRunsRun2            &&
                                       kTPCTempRun2              &&
                                       kTriggerTestsRun2         &&
                                       kFailedRunsRun2           &&
                                       kUnspecifiedTestsRun2     &&
                                       kServerIssuesRun2         &&
                                       kTaggedByRunCoRun2        &&
                                       kPMTTestsRun2             &&
                                       kCRTTestsRun2             &&
                                       kAbsentFromRunCoSheetRun2 &&
                                       kDAQTestsRun2             &&
                                       kTPCTestsRun2             &&
                                       kWireBiasRun2             &&
                                       kShortRunsRun2            &&
                                       kMissingTPCRun2           &&
                                       kInvalidTimingRun2        &&
                                       kAbsentDBRun2             &&
                                       kTestConfigRun2           &&
                                       kMissingPMTRun2           ;

  // Pandora Clear Cosmics per Event
  const ana::SpillCut kDQClearCosmicRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_ClearCosmics_tagged_run2(spill->hdr.run);
    }
  );

  // Number of PMT Flashes per Event
  const ana::SpillCut kDQFlashesRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_Flashes_tagged_run2(spill->hdr.run);
    }
  );

  // Number of PMT Flashes per CC per Event
  const ana::SpillCut kDQFlashesPerCCRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_FlashesPerCC_tagged_run2(spill->hdr.run);
    }
  );

  // Beam-like slices per Event
  const ana::SpillCut kDQBeamLikeRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_dqBeamLike_tagged_run2(spill->hdr.run);
    }
  );

  // Slices per Event
  const ana::SpillCut kDQSlicesRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_sqSlices_tagged_run2(spill->hdr.run);
    }
  );

  // Tracked Hits (any wire plane) per Event
  const ana::SpillCut kDQTrackedHitsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_TrackedHits_tagged_run2(spill->hdr.run);
    }
  );

  // Hits per Track (any wire plane) per Event
  const ana::SpillCut kDQHitsPerTrackRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_HitsPerTrack_tagged_run2(spill->hdr.run);
    }
  );

  // Voltage
  const ana::SpillCut kDQVoltageRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return not is_Voltage_tagged_run2(spill->hdr.run);
    }
  );

  // Quality Metrics
  const ana::SpillCut kDQMetricsRun2 = kDQClearCosmicRun2  &&
                                       kDQFlashesRun2      &&
                                       kDQFlashesPerCCRun2 &&
                                       kDQBeamLikeRun2     &&
                                       kDQSlicesRun2       &&
                                       kDQTrackedHitsRun2  &&
                                       kDQHitsPerTrackRun2 &&
                                       kDQVoltageRun2      ;

  // Good Runs 
  const ana::SpillCut kGoodRunsRun2
  (
    [](const caf::SRSpillProxy* spill)
    {
      return is_goodrun_run2(spill->hdr.run);
    }
  );
} // end ana::icarus namespace
