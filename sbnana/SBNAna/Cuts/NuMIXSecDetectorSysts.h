//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/CAFAna/Core/ISyst.h"

namespace ana
{

  class NuMIXSecDetectorSysts: public ISyst
  {
  public:

    enum DetSystType {
      kFrontIndPlaneGain=0,
      kFrontIndPlaneNoise=1,
      kFrontIndPlaneSignalShape=2,
      kMiddleIndPlaneTransparency=3,
      kSCE=4,
    };

    NuMIXSecDetectorSysts(DetSystType detsyst_type, const std::string& name, const std::string& latexName);

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

  private:

    DetSystType kDetSystType;

  };

  extern const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneGainSyst;
  extern const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneNoiseSyst;
  extern const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneSignalShapeSyst;
  extern const NuMIXSecDetectorSysts kNuMIXSecMiddleIndPlaneTransparencySyst;

}
