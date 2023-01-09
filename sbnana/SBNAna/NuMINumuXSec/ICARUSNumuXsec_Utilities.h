#pragma once

#include "cetlib/search_path.h"

#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TVector3.h"
#include "TMatrixD.h"
#include "TProfile.h"
#include "TFile.h"
#include <iostream>

using namespace std;

namespace ICARUSNumuXsec{

  //==== TODO Use spill info in the future..
  enum GateType
  {
    UNKNOWN = 0,
    BNB = +1,
    BNBOffBeam = -1,
    NUMI = +2,
    NUMIOffBeam = -2,
  };

  class FiducialVolumeTool{

  public:

    FiducialVolumeTool()
    {
      fvCryo0 = {-358.489, -61.94,  // x
                 -181.86, +134.96,  // y
                 -894.95, +894.85};
      fvCryo1 = { +61.94, +358.489,  // x
                 -181.86, +134.96,  // y
                 -894.95, +894.85};
    }

    FidVol fvCryo0, fvCryo1;
    mutable double XMargin, YMargin, ZMarginUp, ZMarginDown;

    bool isContained(double x, double y, double z) const;
    int containedCryo(double x, double y, double z) const;

    int TPCIndex(double x, double y, double z) const;


  }; // END Class FiducialVolumeTool

  class VertexContained : public FiducialVolumeTool{

  public:

    VertexContained(){
      XMargin = 25.;
      YMargin = 25.;
      ZMarginUp = 30.;
      ZMarginDown = 50.;
    }

    static VertexContained& Instance();

  }; // END VertexContained class def

  class TrackContained : public FiducialVolumeTool{

  public:

    TrackContained(){
      XMargin = 10.;
      YMargin = 10.;
      ZMarginUp = 10.;
      ZMarginDown = 10.;
    }

    static TrackContained& Instance();

  }; // END TrackContained class def

  //==== For a given truth particle, find the reco object
  int GetMatchedRecoTrackIndex(const caf::SRSliceProxy* slc, int truth_idx);
  int GetMatchedRecoShowerIndex(const caf::SRSliceProxy* slc, int truth_idx);
  int GetMatchedRecoStubIndex(const caf::SRSliceProxy* slc, int truth_idx);

  //==== Stub calib
  double GetEnergyFromStubCharge(double q);

  class NuMICoordinateTool{

  public:

    NuMICoordinateTool();

    TVector3 GetICARUSCoord(double x, double y, double z) const;

    static NuMICoordinateTool& Instance();

  protected:

    TVector3 NuDirection_NuMI;
    TMatrixD rotMatNtoI = TMatrixD(3,3);
    TMatrixD tranVecNtoI = TMatrixD(3,1);

  }; // END Class NuMICoordinateTool

  class dEdXTemplateTool{

  public:

    dEdXTemplateTool();

    double GetdEdX(double rr, int ptlType) const;
    double GetdEdXErr(double rr, int ptlType) const;

    static dEdXTemplateTool& Instance();

    std::string fTemplateFile;
    std::string fROOTfile;

    TProfile *dedx_range_pro;   ///< proton template
    TProfile *dedx_range_ka;    ///< kaon template
    TProfile *dedx_range_pi;    ///< pion template
    TProfile *dedx_range_mu;    ///< muon template

  };

  class SterileNuTool{

  public:

    SterileNuTool();

    double GetOscProb(double LoE, int i, int f) const;

    static SterileNuTool& Instance();

    double sin2th;
    double m2;

  };

  class CRTPMTMatchingTool{

    public:

    CRTPMTMatchingTool();

    static CRTPMTMatchingTool& Instance();

    void SetGateType(GateType gt) const;
    void SetInTimeRange(double t_min, double t_max) const;
    bool IsInTime(double t_gate) const;
    int GetMatchedCRTHitIndex(
      double opt,
      const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
      int mode
    ) const;
    int GetMatchID(
      double opt,
      const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits
    ) const;

    bool IsNegativeTOF(double timediff) const;

    mutable GateType GT;
    mutable double timecut_min, timecut_max;

  };

}
