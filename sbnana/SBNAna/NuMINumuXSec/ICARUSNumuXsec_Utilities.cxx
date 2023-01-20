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
      auto const& trk = slc->reco.pfp.at(i).trk;
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
      auto const& shw = slc->reco.pfp.at(i).shw;
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
      auto const& stub = slc->reco.stub.at(i);
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

int CRTPMTMatchingTool::GetMatchedCRTHitIndex(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
  int mode) const{

  double mindiff = std::numeric_limits<double>::max();
  int ret=-1;
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& hit = crt_hits.at(i);
    if(hit.plane>=30 && hit.plane<=34){
      double crtt = UseTS0 ? hit.t0 : hit.t1;
      double this_diff = crtt-opt;
      if(fabs(this_diff)<fabs(mindiff)){
        mindiff = this_diff;
        ret = i;
      }
    }
  }

  return ret;

}

int CRTPMTMatchingTool::GetMatchID(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits) const{

  int hasCRTHit = 0;
  static double interval = 0.1;
  int topen = 0, topex = 0, sideen = 0, sideex = 0;
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& hit = crt_hits.at(i);
    double crtt = UseTS0 ? hit.t0 : hit.t1;
    double tof = crtt - opt;

    if(tof<0 && abs(tof)<interval){
      if(hit.plane > 36){
        sideen++;
      }
      else{
        topen++;
      }
    }
    else if(tof>=0 && abs(tof)<interval){
      if(hit.plane > 36){
        sideex++;
      }
      else{
        topex++;
      }
    }

  }

  // hasCRTHit = 0, no matched CRT
  // hasCRTHit = 1, 1 entering from Top CRT
  // hasCRTHit = 2, 1 entering from Side CRT
  // hasCRTHit = 3, 1 entering from Top and exiting to Side CRT
  // hasCRTHit = 4, No entering; 1 exiting to top
  // hasCRTHit = 5, No entering, 1 exiting to side
  // hasCRTHit = 6, Multiple entering
  // hasCRTHit = 7, Multiple entering and exiting to side
  // hasCRTHit = 8, all other cases

  if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0)
    hasCRTHit = 0;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0)
    hasCRTHit = 1;
  else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0)
    hasCRTHit = 2;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1)
    hasCRTHit = 3;
  else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0)
    hasCRTHit = 4;
  else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1)
    hasCRTHit = 5;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0)
    hasCRTHit = 6;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1)
    hasCRTHit = 7;
  else
    hasCRTHit = 8;

  return hasCRTHit;

}

bool CRTPMTMatchingTool::IsNegativeTOF(double timediff) const{
  return (timediff>-0.1 && timediff<0.);
}

