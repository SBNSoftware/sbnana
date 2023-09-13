// B. Howard - 2023
//   howard <at> fnal.gov
// The cuts borrow heavily from NOvA, thanks NOvA!

#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/MultiVar.h"

#include "sbnanaobj/StandardRecord/SRNuMIInfo.h"

namespace ana
{
  double HornCurrentVal ( const caf::SRNuMIInfoProxy& info );
  bool HornCurrentCut ( const caf::SRNuMIInfoProxy& info, const bool stringentMode );

  double POTInSpillVal ( const caf::SRNuMIInfoProxy& info );
  bool POTInSpillCut ( const caf::SRNuMIInfoProxy& info, const bool stringentMode );

  double extrapToLoc ( const double var1, const double loc1, const double var2, const double loc2, const double loc3 );
  std::pair<double,double> BeamPositionAtTargetVal( const caf::SRNuMIInfoProxy& info, const unsigned int runNumber );
  bool BeamPositionAtTargetCut ( const caf::SRNuMIInfoProxy& info, const unsigned int runNumber, const bool stringentMode );

  std::pair<double,double> BeamWidthVal ( const caf::SRNuMIInfoProxy& info );
  bool BeamWidthCut ( const caf::SRNuMIInfoProxy& info );

  // For all spills -- only should return info for first event in subrun
  extern const SpillMultiVar kHornCurrentAll;
  extern const SpillMultiVar kPOTInSpillAll;
  extern const SpillMultiVar kBeamPosHAll;
  extern const SpillMultiVar kBeamPosVAll;
  extern const SpillMultiVar kBeamWidthHAll;
  extern const SpillMultiVar kBeamWidthVAll;

  // For our "triggering" spills
  extern const SpillVar kHornCurrent;
  extern const SpillVar kPOTInSpill;
  extern const SpillVar kBeamPosH;
  extern const SpillVar kBeamPosV;
  extern const SpillVar kBeamWidthH;
  extern const SpillVar kBeamWidthV;

  // Exposure accounting with and without cuts...
  extern const SpillVar kDummyVarForPOTCounting;
  extern const SpillVar kSummedPOT_NuMI_All;
  extern const SpillVar kSummedPOT_NuMI_TRTGTD_All;
  extern const SpillVar kSummedPOT_NuMI_Cuts;
}