#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

#include <algorithm>
#include <functional>

/**
 * ICARUSDataQualityCuts.h
 * Store list of DAQ runs which have known features
 * Runs known at compile time, so list them as const
 * author: hhausner@fnal.gov
 **/

namespace ana::icarus
{
  // helper template lambda to make defining the quality checks easier
  // there is a small overhead defining the functions this way
  // but there are a lot of these, so bare with me
  template <size_t N>
    auto inRunArray (const std::array<unsigned int, N>& runSet)
    {
      return [&runSet](const unsigned int& run) -> bool
      {
        if (std::find(runSet.begin(), runSet.end(), run) != runSet.end())
        {
          return true;
        }
        return false;
      };
    }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Prefiltered Runs: https://docs.google.com/spreadsheets/d/1Kra6eIflTKS_sMghBqgpy1h86Z8WKkibZLrLnDAdWaQ
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Laser Runs
  constexpr std::array<unsigned int, 13>
    laserRuns_run2{ 9304,  9305,  9306,  9628,  9629,  9630,  9656,  9657,  9866,  9936,
                    9937,  9938,  9939                                                 };
  auto isLaserRun_run2 = inRunArray(laserRuns_run2);

  // TPC Temp
  constexpr std::array<unsigned int, 7> 
    tpcTempIssues_run2{ 9320,  9321,  9322,  9323,  9324,  9325,  9326};
  auto hasTPCTempIssue_run2 = inRunArray(tpcTempIssues_run2);

  // Trigger Tests
  constexpr std::array<unsigned int, 168>
    triggerTests_run2{ 9348,  9349,  9350,  9351,  9352,  9527,  9528,  9635,  9636,  9637,
                       9638,  9639,  9640,  9641,  9650,  9651,  9652,  9653,  9654,  9655,
                       9659,  9660,  9661,  9662,  9663,  9664,  9665,  9666,  9667,  9668,
                       9669,  9670,  9671,  9676,  9677,  9678,  9679,  9680,  9681,  9682,
                       9683,  9684,  9685,  9686,  9736,  9737,  9738,  9739,  9740,  9741,
                       9767,  9768,  9769,  9770,  9771,  9772,  9773,  9774,  9775,  9776,
                       9777,  9778,  9787,  9789,  9790,  9798,  9799,  9800,  9801,  9802,
                       9803,  9804,  9805,  9806,  9872,  9873,  9874,  9875,  9876,  9877,
                       9878,  9879,  9880,  9881,  9882,  9883,  9884,  9885,  9886,  9887,
                       9888,  9889,  9890,  9891,  9898,  9899,  9900,  9901,  9902,  9903,
                       9904,  9905,  9906,  9907,  9908,  9909,  9910,  9957,  9958,  9996,
                       9997,  9998,  9999, 10000, 10001, 10002, 10003, 10004, 10005, 10006,
                      10007, 10008, 10016, 10017, 10018, 10019, 10020, 10021, 10022, 10023,
                      10024, 10025, 10026, 10027, 10028, 10029, 10035, 10036, 10037, 10038,
                      10039, 10041, 10042, 10043, 10044, 10045, 10046, 10047, 10049, 10050,
                      10051, 10052, 10053, 10069, 10070, 10071, 10072, 10073, 10074, 10075,
                      10076, 10077, 10078, 10079, 10080, 10081, 10082, 10083              };
  auto isTriggerTest_run2 = inRunArray(triggerTests_run2);

  // Failed Runs
  constexpr std::array<unsigned int, 5>
    failedRuns_run2{ 9381,  9382,  9446,  9447,  9561};
  auto isFailedRun_run2 = inRunArray(failedRuns_run2);

  // Unspecified Tests
  constexpr std::array<unsigned int, 55>
    unspecifiedTests_run2{ 9395,  9396,  9397,  9398,  9399,  9400,  9401,  9402,  9403,  9404,
                           9405,  9406,  9407,  9408,  9410,  9411,  9416,  9434,  9440,  9459,
                           9465,  9466,  9467,  9476,  9503,  9530,  9567,  9808,  9809,  9810,
                           9811,  9812,  9813,  9815,  9816,  9817,  9818,  9819,  9820,  9821,
                           9822,  9823,  9824,  9825,  9826,  9827,  9828,  9829,  9830,  9831,
                           9832,  9833,  9864,  9865,  9895                                   };
  auto isUnspecifiedTest_run2 = inRunArray(unspecifiedTests_run2);

  // Server Issues
  constexpr std::array<unsigned int, 2>
    serverIssues_run2{ 9413,  9414};
  auto hasServerIssues_run2 = inRunArray(serverIssues_run2);


  // Tagged By RunCo
  constexpr std::array<unsigned int, 12>
    taggedByRunCo_run2{ 9433,  9449,  9500,  9501,  9502,  9596,  9600,  9601,  9759, 9952,
                        9969,  9976                                                       };
  auto isTaggedByRunCo = inRunArray(taggedByRunCo_run2); 

  // PMT Tests
  constexpr std::array<unsigned int, 36>
    pmtTests_run2{ 9451,  9452,  9453,  9463,  9468,  9469,  9470,  9471,  9479,  9480,
                   9492,  9493,  9494,  9505,  9506,  9507,  9508,  9509,  9510,  9511,
                   9512,  9529,  9572,  9573,  9574,  9575,  9576,  9581,  9604,  9605,
                   9606,  9607,  9608,  9609,  9845,  9846                            };
  auto isPMTTest_run2 = inRunArray(pmtTests_run2);

  // CRT Tests
  constexpr std::array<unsigned int, 20>
    crtTests_run2{ 9454,  9455,  9456,  9462,  9496,  9497,  9498,  9519,  9520,  9521,
                   9524,  9525,  9526,  9579,  9603,  9912,  9992,  9993,  9994,  9995};
  auto isCRTTest_run2 = inRunArray(crtTests_run2);

  // Absent from RunCo Sheet
  constexpr std::array<unsigned int, 35>
    notInSheet_run2{ 9461,  9483,  9484,  9485,  9486,  9487,  9488,  9489,  9490,  9491,
                     9535,  9536,  9537,  9538,  9539,  9540,  9541,  9542,  9543,  9544,
                     9545,  9546,  9547,  9548,  9549,  9550,  9551,  9552,  9553,  9554,
                     9555,  9556,  9557,  9577,  9578                                   };
  auto isNotInSheet_run2 = inRunArray(notInSheet_run2);

  // DAQ Tests
  constexpr std::array<unsigned int, 14>
    daqTests_run2{ 9522,  9523,  9632,  9633,  9643,  9644,  9645,  9706,  9707,  9708,
                   9709,  9856,  9857,  9930                                          };
  auto isDAQTest_run2 = inRunArray(daqTests_run2);

  // TPC Tests
  constexpr std::array<unsigned int, 21>
    tpcTests_run2{ 9591,  9592,  9701,  9702,  9710,  9711,  9712,  9713,  9718,  9719,
                   9927,  9928,  9990, 10011, 10012, 10013, 10014, 10031, 10032, 10033,
                  10034                                                               };
  auto isTPCTest_run2 = inRunArray(tpcTests_run2);

  // Power Issues at FD
  constexpr std::array<unsigned int, 13>
    powerIssues_run2{ 9611,  9612,  9613,  9614,  9615,  9616,  9617,  9618,  9619,  9620,
                      9621,  9622,  9623                                                 };
  auto hasPowerIssue_run2 = inRunArray(powerIssues_run2);

  // WireBias PS Problems
  constexpr std::array<unsigned int, 5>
    wireBiasIssues_run2{ 9963,  9964,  9965,  9966,  9967};
  auto hasWireBiasIssue_run2 = inRunArray(wireBiasIssues_run2);

  // Short Runs
  constexpr std::array<unsigned int, 37>
    shortRuns_run2{ 9315,  9319,  9428,  9429,  9430,  9432,  9442,  9443,  9444,  9457,
                    9464,  9475,  9495,  9513,  9514,  9531,  9532,  9559,  9566,  9571,
                    9584,  9585,  9625,  9634,  9687,  9722,  9727,  9782,  9784,  9785,
                    9786,  9836,  9839,  9842,  9931,  9933,  9934                     };
  auto isShortRun_run2 = inRunArray(shortRuns_run2);

  // Missing TPC Components
  constexpr std::array<unsigned int, 21>
    missingTPC_run2{ 9331,  9355,  9427,  9515,  9516,  9517,  9586,  9624,  9673,  9674,
                     9720,  9748,  9749,  9751,  9754,  9756,  9757,  9760,  9761,  9766,
                     9867                                                               };
  auto isMissingTPC_run2 = inRunArray(missingTPC_run2);

  // Invalid Timing Information
  constexpr std::array<unsigned int, 2>
    missingTimingInfo_run2{ 9334,  9797};
  auto isMissingTimingInfo_run2 = inRunArray(missingTimingInfo_run2);

  // Absent from Database
  constexpr std::array<unsigned int, 22>
    notInDatabase_run2{ 9368,  9369,  9370,  9371,  9372,  9373,  9374,  9375,  9376,  9377,
                        9378,  9379,  9779,  9780,  9781,  9852,  9853,  9913,  9962,  9973,
                        9975,  9988                                                        };
  auto isNotInDatabase_run2 = inRunArray(notInDatabase_run2);


  // Test Configuration (Test_thr390_LockTemp_true_00001 only)
  constexpr std::array<unsigned int, 10>
    testConfig_run2{ 9417,  9418,  9419,  9420,  9421,  9422,  9423,  9424,  9425,  9426};
  auto isTestConfig_run2 = inRunArray(testConfig_run2);

  // Missing PMT Components
  constexpr std::array<unsigned int, 2>
    missingPMT_run2{ 9431,  9814};
  auto isMissingPMT_run2 = inRunArray(missingPMT_run2);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data quality metrics: https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=40006
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Pandora Clear Cosmics per Event
  constexpr std::array<unsigned int, 1>
    dqClearCosmics_run2{ 9984 };
  auto is_ClearCosmics_tagged_run2 = inRunArray(dqClearCosmics_run2);

  // Number of PMT Flashes per Event
  constexpr std::array<unsigned int, 10>
    dqFlashes_run2{ 9336,  9533,  9534,  9558,  9560,  9562,  9564,  9568,  9580,  9582};
  auto is_Flashes_tagged_run2 = inRunArray(dqFlashes_run2);

  // Number of PMT Flashes per CC per Event
  constexpr std::array<unsigned int, 9>
    dqFlashesPerCC_run2{ 9533,  9534,  9558,  9560,  9562,  9564,  9568,  9580,  9582};
  auto is_FlashesPerCC_tagged_run2 = inRunArray(dqFlashesPerCC_run2);

  // Beam-like slices per Event
  constexpr std::array<unsigned int, 3>
    dqBeamLike_run2{ 9345,  9642,  9972};
  auto is_dqBeamLike_tagged_run2 = inRunArray(dqBeamLike_run2);

  // Slices per Event
  constexpr std::array<unsigned int, 1>
    dqSlices_run2{10060};
  auto is_sqSlices_tagged_run2 = inRunArray(dqSlices_run2);

  // Tracked Hits (any wire plane) per Event
  constexpr std::array<unsigned int, 3>
    dqTrackedHits_run2{ 9968,  9984, 10060};
  auto is_TrackedHits_tagged_run2 = inRunArray(dqTrackedHits_run2);

  // Hits per Track (any wire plane) per Event
  constexpr std::array<unsigned int, 4>
    dqHitsPerTrack_run2{ 9580,  9968,  9983,  9984};
  auto is_HitsPerTrack_tagged_run2 = inRunArray(dqHitsPerTrack_run2);

  // Voltage
  constexpr std::array<unsigned int, 3>
    dqVoltage_run2{ 9610,  9983,  9984};
  auto is_Voltage_tagged_run2 = inRunArray(dqVoltage_run2);
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Good Runs
  // This almost (but not exactly) all runs which are not prefiltered or tagged by the above
  // Run2 Exceptions: 
  //   - 9435 is excluded as it is the handscan run
  //   - 9610 is included as we believe the power glitch only effected the end of the run
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  constexpr std::array<unsigned int, 230>
    goodruns_run2{ 9301,  9302,  9303,  9307,  9308,  9309,  9310,  9311,  9312,  9313,
                   9314,  9316,  9317,  9318,  9327,  9328,  9329,  9330,  9332,  9333,
                   9335,  9337,  9338,  9339,  9340,  9341,  9342,  9343,  9344,  9346,
                   9347,  9353,  9354,  9356,  9357,  9358,  9359,  9360,  9361,  9362,
                   9363,  9364,  9365,  9366,  9367,  9380,  9383,  9384,  9385,  9386,
                   9387,  9388,  9389,  9390,  9391,  9392,  9393,  9394,  9409,  9412,
                   9415,  9435,  9436,  9437,  9438,  9439,  9441,  9445,  9448,  9450,
                   9458,  9460,  9472,  9473,  9474,  9477,  9478,  9481,  9482,  9499,
                   9504,  9518,  9563,  9565,  9569,  9570,  9583,  9587,  9588,  9589,
                   9590,  9593,  9594,  9595,  9597,  9598,  9599,  9602,  9610,  9626,
                   9627,  9631,  9647,  9648,  9649,  9658,  9672,  9675,  9688,  9689,
                   9690,  9691,  9692,  9693,  9694,  9695,  9696,  9697,  9698,  9699,
                   9700,  9703,  9704,  9705,  9714,  9715,  9716,  9717,  9721,  9723,
                   9724,  9725,  9726,  9728,  9729,  9730,  9731,  9732,  9733,  9734,
                   9735,  9743,  9744,  9745,  9746,  9747,  9750,  9752,  9753,  9755,
                   9758,  9762,  9763,  9764,  9765,  9783,  9788,  9791,  9792,  9793,
                   9794,  9795,  9796,  9807,  9834,  9835,  9837,  9838,  9840,  9841,
                   9844,  9847,  9849,  9851,  9854,  9855,  9860,  9862,  9868,  9870,
                   9892,  9894,  9896,  9897,  9914,  9917,  9919,  9921,  9922,  9924,
                   9925,  9926,  9929,  9932,  9935,  9940,  9941,  9942,  9944,  9945,
                   9946,  9949,  9950,  9951,  9953,  9954,  9956,  9959,  9960,  9961,
                   9970,  9971,  9974,  9977,  9979,  9981,  9982,  9986, 10054, 10059,
                  10061, 10062, 10064, 10065, 10066, 10067, 10084, 10085, 10096, 10097};
  auto is_goodrun_run2 = inRunArray(goodruns_run2);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // CAF cuts
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // does the run pass our prefiltering?
  extern const ana::SpillCut kLaserRunsRun2;
  extern const ana::SpillCut kTPCTempRun2;
  extern const ana::SpillCut kTriggerTestsRun2;
  extern const ana::SpillCut kFailedRunsRun2;
  extern const ana::SpillCut kUnspecifiedTestsRun2;
  extern const ana::SpillCut kServerIssuesRun2;
  extern const ana::SpillCut kTaggedByRunCoRun2;
  extern const ana::SpillCut kPMTTestsRun2;
  extern const ana::SpillCut kCRTTestsRun2;
  extern const ana::SpillCut kAbsentFromRunCoSheetRun2;
  extern const ana::SpillCut kDAQTestsRun2;
  extern const ana::SpillCut kTPCTestsRun2;
  extern const ana::SpillCut kWireBiasRun2;
  extern const ana::SpillCut kShortRunsRun2;
  extern const ana::SpillCut kMissingTPCRun2;
  extern const ana::SpillCut kInvalidTimingRun2;
  extern const ana::SpillCut kAbsentDBRun2;
  extern const ana::SpillCut kTestConfigRun2;
  extern const ana::SpillCut kMissingPMTRun2;
  extern const ana::SpillCut kPrefilterRun2;

  // does the run pass our quality metrics?
  extern const ana::SpillCut kDQClearCosmicRun2;
  extern const ana::SpillCut kDQFlashesRun2;
  extern const ana::SpillCut kDQFlashesPerCCRun2;
  extern const ana::SpillCut kDQBeamLikeRun2;
  extern const ana::SpillCut kDQSlicesRun2;
  extern const ana::SpillCut kDQTrackedHitsRun2;
  extern const ana::SpillCut kDQHitsPerTrackRun2;
  extern const ana::SpillCut kDQVoltageRun2;
  extern const ana::SpillCut kDQMetricsRun2;

  // is the run in our good runs list?
  extern const ana::SpillCut kGoodRunsRun2;

} // end ana::icarus namespace
