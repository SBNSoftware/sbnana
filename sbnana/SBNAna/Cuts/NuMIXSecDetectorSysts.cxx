#include "sbnana/SBNAna/Cuts/NuMIXSecDetectorSysts.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include <iostream>

namespace ana {

  NuMIXSecDetectorSysts::NuMIXSecDetectorSysts(DetSystType detsyst_type, const std::string& name, const std::string& latexName):
    ISyst(name, latexName),
    kDetSystType(detsyst_type)
  {

  }
  void NuMIXSecDetectorSysts::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {

    int RecoMuonIdx = kNuMIMuonCandidateIdx(sr);
    if(RecoMuonIdx<0) return;
    int RecoProtonIdx = kNuMIProtonCandidateIdx(sr);
    if(RecoProtonIdx<0) return;

    if(kDetSystType==kFrontIndPlaneGain){
      double RecoProtonP = kNuMIProtonCandidateRecoP(sr);
      weight *= 1. + sigma * ( RecoProtonP<0.5 ? 0.10 : 0.05 );
    }
    else if(kDetSystType==kFrontIndPlaneNoise){
      // +1: reduce the rate by 10%
      //  0: CV
      double this_sigma = sigma<0. ? 0. : sigma;
      weight *= 1. + this_sigma * (-0.10);
    }
    else if(kDetSystType==kFrontIndPlaneSignalShape){
      // +1: reduce the rate by 10%
      //  0: CV
      double this_sigma = sigma<0. ? 0. : sigma;
      weight *= 1. + this_sigma * (-0.10);
    }
    else if(kDetSystType==kMiddleIndPlaneTransparency){
      weight *= 1. + sigma * 0.05;
    }
    else{

    }


  }

  void NuMIXSecDetectorSysts::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

  }

  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneGainSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneGain,
    "NuMIXSecFrontIndPlaneGainSyst",
    "Front ind. plane gain #pm10%"
  );
  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneNoiseSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneNoise,
    "NuMIXSecFrontIndPlaneNoiseSyst",
    "Front ind. plane noise +10%"
  );
  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneSignalShapeSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneSignalShape,
    "NuMIXSecFrontIndPlaneSignalShapeSyst",
    "Front ind. plane signal shape"
  );
  const NuMIXSecDetectorSysts kNuMIXSecMiddleIndPlaneTransparencySyst(
    NuMIXSecDetectorSysts::kMiddleIndPlaneTransparency,
    "NuMIXSecMiddleIndPlaneTransparencySyst",
    "Middle ind. plane transparency"
  );


} // end namespace ana
