#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"

using namespace ICARUSNumuXsec;

VertexContained& VertexContained::Instance(){
  static VertexContained vc;
  return vc;
}

TrackContained& TrackContained::Instance(){
  static TrackContained tc;
  return tc;
}

bool FiducialVolumeTool::isContained(double x, double y, double z) const {

  int _containedCryo = containedCryo(x,y,z);
  return (_containedCryo==0 || _containedCryo==1);

}

int FiducialVolumeTool::containedCryo(double x, double y, double z) const {

  int out=-1;
  //==== Cryo0
  if( x <= 0 ){
    if ( !isnan(x) &&
         ( x < fvCryo0.xmax - XMargin && x > fvCryo0.xmin + XMargin ) &&
         !isnan(y) &&
         ( y > fvCryo0.ymin + YMargin && y < fvCryo0.ymax - YMargin ) &&
         !isnan(z) &&
         ( z > fvCryo0.zmin + ZMarginUp && z < fvCryo0.zmax - ZMarginDown )
       ) out = 0;
  }
  //==== Cryo1
  else{
    if ( !isnan(x) &&
         ( x > fvCryo1.xmin + XMargin && x < fvCryo1.xmax - XMargin ) &&
         !isnan(y) &&
         ( y > fvCryo0.ymin + YMargin && y < fvCryo0.ymax - YMargin ) &&
         !isnan(z) &&
         ( z > fvCryo0.zmin + ZMarginUp && z < fvCryo0.zmax - ZMarginDown )
       ) out = 1;
  }
  return out;

}

int FiducialVolumeTool::TPCIndex(double x, double y, double z) const {

  double x_cath = x<=0 ? (fvCryo0.xmax+fvCryo0.xmin)/2. : (fvCryo1.xmax+fvCryo1.xmin)/2.;
  double x_from_cath = x-x_cath;
  if(x_from_cath<0 && z<0) return 0;
  else if(x_from_cath<0 && z>0) return 1;
  else if(x_from_cath>0 && z<0) return 2;
  else if(x_from_cath>0 && z>0) return 3;
  else return 4;

}

//==== For a given truth particle, find the reco object
//====   Track : return the longest matched track
int ICARUSNumuXsec::GetMatchedRecoTrackIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    double LMax(-999.);
    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){
      const auto& pfp = slc->reco.pfp.at(i);
      if(pfp.trackScore<0.5) continue;
      const auto& trk = pfp.trk;
      if( trk.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && trk.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        if(trk.len > LMax){
          PTrackInd = i;
          LMax = trk.len;
        }
      }
    }
    return PTrackInd;
  }
  else{
    return -1;
  }

}

//====   Shower : return (TODO energetic?) shower
int ICARUSNumuXsec::GetMatchedRecoShowerIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){
      const auto& pfp = slc->reco.pfp.at(i);
      if(pfp.trackScore>=0.5) continue;
      const auto& shw = pfp.shw;
      if( shw.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && shw.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        PTrackInd = i;
      }
    }
    return PTrackInd;
  }
  else{
    return -1;
  }

}

int ICARUSNumuXsec::GetMatchedRecoStubIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    for(std::size_t i(0); i < slc->reco.stub.size(); ++i){
      const auto& stub = slc->reco.stub.at(i);
      if( stub.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && stub.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        PTrackInd = i;
      }
    }
    return PTrackInd;
  }
  else{
    return -1;
  }

}

double ICARUSNumuXsec::GetEnergyFromStubCharge(double q){
  static double p0 = -0.00236892;
  static double p1 = 0.383538;
  if(q<=0.) return -999.;
  else return (q*23.6e-9-p0)/p1; // to GeV
}

NuMICoordinateTool::NuMICoordinateTool(){

  NuDirection_NuMI.SetXYZ(3.94583e-01, 4.26067e-02, 9.17677e-01);

  TMatrixDRow(rotMatNtoI,0) = {0.921035925, 0.022715103, 0.388814672};
  TMatrixDRow(rotMatNtoI,1) = {0., 0.998297825, -0.058321970};
  TMatrixDRow(rotMatNtoI,2) = {-0.389477631, 0.053716629, 0.919468161};
  rotMatNtoI.Print();

  TMatrixDColumn(tranVecNtoI,0) = {-315.120380, -33.644912, -733.632532};
  tranVecNtoI.Print();

}

TVector3 NuMICoordinateTool::GetICARUSCoord(double x, double y, double z) const {

  TMatrixD coordN(3,1);
  TMatrixDColumn(coordN,0) = {x, y, z};

  TMatrixD ret = (rotMatNtoI*coordN+tranVecNtoI);

  return TVector3(ret(0,0), ret(1,0), ret(2,0));

}

NuMICoordinateTool& NuMICoordinateTool::Instance(){
  static NuMICoordinateTool nct;
  return nct;
}

dEdXTemplateTool::dEdXTemplateTool(){

  std::cout << "[dEdXTemplateTool::dEdXTemplateTool] Setting up.." << std::endl;

  cet::search_path sp("FW_SEARCH_PATH");

  std::string fTemplateFile = "dEdxrestemplates.root";

  sp.find_file(fTemplateFile, fROOTfile);

  TFile *file = TFile::Open(fROOTfile.c_str());
  dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
  dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
  dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
  dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

}

double dEdXTemplateTool::GetdEdX(double rr, int ptlType) const {

  int bin = dedx_range_pro->FindBin(rr);
  if(bin>dedx_range_pro->GetNbinsX()) bin = dedx_range_pro->GetNbinsX()-1;
  if(bin>=1&&bin<=dedx_range_pro->GetNbinsX()){

    double bincpro = dedx_range_pro->GetBinContent(bin);
    if (bincpro<1e-6){//for 0 bin content, using neighboring bins
      bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2;
    }
    double bincka = dedx_range_ka->GetBinContent(bin);
    if (bincka<1e-6){
      bincka = (dedx_range_ka->GetBinContent(bin-1)+dedx_range_ka->GetBinContent(bin+1))/2;
    }
    double bincpi = dedx_range_pi->GetBinContent(bin);
    if (bincpi<1e-6){
      bincpi = (dedx_range_pi->GetBinContent(bin-1)+dedx_range_pi->GetBinContent(bin+1))/2;
    }
    double bincmu = dedx_range_mu->GetBinContent(bin);
    if (bincmu<1e-6){
      bincmu = (dedx_range_mu->GetBinContent(bin-1)+dedx_range_mu->GetBinContent(bin+1))/2;
    }

    if(ptlType==0) return bincpro;
    else if(ptlType==1) return bincka;
    else if(ptlType==2) return bincpi;
    else if(ptlType==3) return bincmu;
    else{
      std::cout << "[dEdXTemplateTool::GetdEdX] Wrong ptlType : " << ptlType << std::endl;
      abort();
      return -1.;
    }

  }
  else{
    std::cout << "[dEdXTemplateTool::GetdEdX] rr = " << rr << ", bin = " << bin << std::endl;
    return -1.;
  }

}

double dEdXTemplateTool::GetdEdXErr(double rr, int ptlType) const {

  int bin = dedx_range_pro->FindBin(rr);
  if(bin>dedx_range_pro->GetNbinsX()) bin = dedx_range_pro->GetNbinsX()-1;
  if(bin>=1&&bin<=dedx_range_pro->GetNbinsX()){

    double binepro = dedx_range_pro->GetBinError(bin);
    if (binepro<1e-6){
      binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2;
    }
    double bineka = dedx_range_ka->GetBinError(bin);
    if (bineka<1e-6){
      bineka = (dedx_range_ka->GetBinError(bin-1)+dedx_range_ka->GetBinError(bin+1))/2;
    }
    double binepi = dedx_range_pi->GetBinError(bin);
    if (binepi<1e-6){
      binepi = (dedx_range_pi->GetBinError(bin-1)+dedx_range_pi->GetBinError(bin+1))/2;
    }
    double binemu = dedx_range_mu->GetBinError(bin);
    if (binemu<1e-6){
      binemu = (dedx_range_mu->GetBinError(bin-1)+dedx_range_mu->GetBinError(bin+1))/2;
    }

    if(ptlType==0) return binepro;
    else if(ptlType==1) return bineka;
    else if(ptlType==2) return binepi;
    else if(ptlType==3) return binemu;
    else{
      std::cout << "[dEdXTemplateTool::GetdEdXErr] Wrong ptlType : " << ptlType << std::endl;
      abort();
      return -1.;
    }

  }
  else{
    std::cout << "[dEdXTemplateTool::GetdEdX] rr = " << rr << ", bin = " << bin << std::endl;
    return -1.;
  }

}

dEdXTemplateTool& dEdXTemplateTool::Instance(){
  static dEdXTemplateTool dedxtt;
  return dedxtt;
}

SterileNuTool::SterileNuTool(){

  sin2th = 0.10;
  m2 = 7.3;

}

double SterileNuTool::GetOscProb(double LoE, int i, int f) const{

  //==== flav : 0/1/2 = e/m/t

  return 1. - (sin2th*sin2th) * pow( TMath::Sin(1.27 * m2 * LoE), 2 );
  
}

SterileNuTool& SterileNuTool::Instance(){
  static SterileNuTool snt;
  return snt;
}

CRTPMTMatchingTool::CRTPMTMatchingTool(){
  std::cout << "[CRTPMTMatchingTool::CRTPMTMatchingTool] called" << std::endl;
  UseTS0 = false;
  Debug = false;
}

CRTPMTMatchingTool& CRTPMTMatchingTool::Instance(){
  static CRTPMTMatchingTool cpmt;
  return cpmt;
}

void CRTPMTMatchingTool::SetGateType(GateType gt) const {
  GT = gt;
  if(abs(GT)==1){
    timecut_min = 0.;
    timecut_max = 2.2;
  }
  else if(abs(GT)==2){
    timecut_min = 0.;
    timecut_max = 10.1;
  }
  else{
    std::cout << "[CRTPMTMatchingTool] Wrong gate type = " << gt << std::endl;
    abort();
  }

  std::cout << "[CRTPMTMatchingTool::SetGateType] GT = " << GT << ", timecut_min = " << timecut_min << ", timecut_max = " << timecut_max << std::endl;

}

void CRTPMTMatchingTool::SetInTimeRange(double t_min, double t_max) const {
  timecut_min = t_min;
  timecut_max = t_max;
  std::cout << "[CRTPMTMatchingTool::SetInTimeRange] timecut_min = " << timecut_min << ", timecut_max = " << timecut_max << std::endl;
}

bool CRTPMTMatchingTool::IsInTime(double t_gate) const{
  //std::cout << "timecut_min = " << timecut_min << ", t_gate = " << t_gate << ", timecut_max = " << timecut_max << std::endl;
  return ( timecut_min<=t_gate && t_gate<=timecut_max );
}

// mode = 0 : Top only
// mode = 1 : Side only
// mode = 2 : Top+Side
std::vector<int> CRTPMTMatchingTool::GetMatchedCRTHitIndex(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
  int mode) const{

  static double interval = 0.1;

  //double mindiff = std::numeric_limits<double>::max();
  std::vector<int> rets = {};
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& hit = crt_hits.at(i);
    bool IsTopCRT = hit.plane>=30 && hit.plane<=39;
    bool IsSideCRT = hit.plane>=40 && hit.plane<=49;
    bool crtSelection = false;
    if(mode==0) crtSelection = IsTopCRT;
    if(mode==1) crtSelection = IsSideCRT;
    if(mode==2) crtSelection = IsTopCRT||IsSideCRT;
    if( crtSelection ){
      double crtt = UseTS0 ? hit.t0 : hit.t1;
      double this_diff = crtt-opt;
      if(abs(this_diff)<interval) rets.push_back(i);
/*
      if(fabs(this_diff)<fabs(mindiff)){
        mindiff = this_diff;
        ret = i;
      }
*/
    }
  }

  return rets;

}

int CRTPMTMatchingTool::GetMatchID(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
  int mode) const{

  int hasCRTHit = -1;
  static double interval = 0.1;
  int TopEn = 0, TopEx = 0, SideEn = 0, SideEx = 0;
  int SideSouthEn = 0;
  int SideEastEn = 0;
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& hit = crt_hits.at(i);
    double crtt = UseTS0 ? hit.t0 : hit.t1;
    double tof = crtt - opt;

    if(tof<0 && abs(tof)<interval){
      if(hit.plane >= 40 && hit.plane<=49){
        SideEn++;
        if(hit.plane==46){
          SideSouthEn++;
        }
        if(hit.plane>=43&&hit.plane<=45){
          SideEastEn++;
        }
      }
      else if(hit.plane >= 30 && hit.plane <= 39){
        TopEn++;
      }
      else{
      }
    }
    else if(tof>=0 && abs(tof)<interval){
      if(hit.plane >= 40 && hit.plane<=49){
        SideEx++;
      }
      else if(hit.plane >= 30 && hit.plane <= 39){
        TopEx++;
      }
      else{
      }
    }

  }

  int OtherSideEn = SideEn - SideSouthEn - SideEastEn;

  // Top+Side
  if(mode==0){

    // -- No matching
    // hasCRTHit = -1, no matched CRT
    // -- 0 entering, but somthing exiting cases (e.g., exiting muon from nu)
    // hasCRTHit = 0, 0 entering, something exiting
    // -- Top-entering cases (Cosmic)
    // hasCRTHit = 1, 1 entering from Top CRT, nothing else
    // hasCRTHit = 2, 1 entering from Top CRT, >=1 exiting to side
    // -- South-entering cases (BNB or NuMI DIRT)
    // hasCRTHit = 3, 1 entering from South-Side CRT (no East-Side entering), nothing else
    // hasCRTHit = 4, 1 entering from South-Side CRT (no East-Side entering), exiting to somewhere
    // -- East-entering cases (NuMI DIRT)
    // hasCRTHit = 5, 1 entering from East-Side CRT (no South-Side entering), nothing else
    // hasCRTHit = 6, 1 entering from Esat-Side CRT (no South-Side entering), exiting to somewhere
    // -- Other side wall entering cases (Cosmic)
    // hasCRTHit = 7, 1 entering from Other-Side CRT, nothing else
    // hasCRTHit = 8, 1 entering from Other-Side CRT, exiting to somewhere
    // -- Multiple top-entering
    // hasCRTHit = 9, >=1 entering from Top CRT, nothing else
    // hasCRTHit = 10, >=1 entering from Top CRT, >=1 exiting to side
    // --- Top-Side entering
    // hasCRTHit = 11, >=1 entering from Top CRT, >=1 entering from side
    // --- TEST
    // hasCRTHit = 12, >=1 entering top and >=1 exiting top

    if(TopEn==0 && SideEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = -1;

    else if(TopEn==0 && SideEn==0 && TopEx+SideEx>=1)
      hasCRTHit = 0;

    else if(TopEn==1 && SideEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 1;
    else if(TopEn==1 && SideEn==0 && TopEx==0 && SideEx>=1)
      hasCRTHit = 2;

    else if(TopEn==0 && SideSouthEn==1 && SideEastEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 3;
    else if(TopEn==0 && SideSouthEn==1 && SideEastEn==0 && TopEx+SideEx>=1)
      hasCRTHit = 4;

    else if(TopEn==0 && SideEastEn==1 && SideSouthEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 5;
    else if(TopEn==0 && SideEastEn==1 && SideSouthEn==0 && TopEx+SideEx>=1)
      hasCRTHit = 6;

    else if(TopEn==0 && OtherSideEn==1 && TopEx==0 && SideEx==0)
      hasCRTHit = 7;
    else if(TopEn==0 && OtherSideEn==1 && TopEx+SideEx>=1)
      hasCRTHit = 8;

    else if(TopEn>=1 && SideEn==0 && TopEx==0 && SideEx==0)
      hasCRTHit = 9;
    else if(TopEn>=1 && SideEn==0 && TopEx==0 && SideEx>=1)
      hasCRTHit = 10;

    else if(TopEn>=1 && SideEn>=1)
      hasCRTHit = 11;

    else if(TopEn>=1 && TopEx>=1)
      hasCRTHit = 12;

    else{
      if(Debug) printf("(TopEn, TopEx, SideEn, SideEx) = (%d, %d, %d, %d)\n",TopEn, TopEx, SideEn, SideEx);
      hasCRTHit = 13;
    }

  }
  // Top only
  else if(mode==1){

    // hasCRTHit = 0, no matched CRT
    // hasCRTHit = 1, One top entering (no Top-exiting)
    // hasCRTHit = 2, One top exiting (no Top-entering)
    // hasCRTHit = 3, Multiple top entering (no top-exiting)
    // hasCRTHit = 4, Multiple top exiting (no top-entering)
    // hasCRTHit = 5, All other cases; entering>=1 && exiting>=1

    if(TopEn==0 && TopEx==0)
      hasCRTHit = 0;
    else if(TopEn==1 && TopEx==0)
      hasCRTHit = 1;
    else if(TopEn==0 && TopEx==1)
      hasCRTHit = 2;
    else if(TopEn>=1 && TopEx==0)
      hasCRTHit = 3;
    else if(TopEn==0 && TopEx>=1)
      hasCRTHit = 4;
    else
      hasCRTHit = 5;

  }
  // side only
  else if(mode==2){

    // hasCRTHit = 0, no matched CRT
    // hasCRTHit = 1, One side entering (no side-exiting)
    // hasCRTHit = 2, One side exiting (no side-entering)
    // hasCRTHit = 3, Multiple side entering (no side-exiting)
    // hasCRTHit = 4, Multiple side exiting (no side-entering)
    // hasCRTHit = 5, All other cases; entering>=1 && exiting>=1

    if(SideEn==0 && SideEx==0)
      hasCRTHit = 0;
    else if(SideEn==1 && SideEx==0)
      hasCRTHit = 1;
    else if(SideEn==0 && SideEx==1)
      hasCRTHit = 2;
    else if(SideEn>=1 && SideEx==0)
      hasCRTHit = 3;
    else if(SideEn==0 && SideEx>=1)
      hasCRTHit = 4;
    else
      hasCRTHit = 5;

  }
  else{

  }

  return hasCRTHit;

}

bool CRTPMTMatchingTool::IsNegativeTOF(double timediff) const{
  return (timediff>-0.1 && timediff<0.);
}

