#include "sbnana/SBNAna/Cuts/NuMIXSecSysts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecTruthVars.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include <iostream>
#include "TMath.h"
#include "TFile.h"

namespace ana {

  //---------------------------------------------------------------------
  // Single-pion production CV correction and systematics

  bool IsSPP(const caf::SRTrueInteractionProxy *sr){

/*
    // Below using NUANCE code
    const auto& NUANCECode = sr->genie_inttype;
    bool IsCCSinglePiPlus = false;
    if(NUANCECode == caf::genie_interaction_type_::kResCCNuProtonPiPlus){
      IsCCSinglePiPlus = true;
    }
    if(NUANCECode == caf::genie_interaction_type_::kResCCNuNeutronPiPlus){
      IsCCSinglePiPlus = true;
    }

    return IsCCSinglePiPlus;
*/

    // Below using event record
    // Check process

    int genie_n_photons = 0;
    int genie_n_mesons = 0;
    int nPip = 0;
    for(const auto& prim: sr->prim){

      // only primary
      if ( prim.start_process != 0 ) continue;

      int pdg = prim.pdg;
      int energy = prim.genE * 1000.; // GeV->MeV

      if (pdg == 22 && energy > 10.0) {
        genie_n_photons++;
      }
      if (abs(pdg) == 211 || //pi+-
               pdg == 111 ||  // pi0
               abs(pdg) == 321 || // K-
               abs(pdg) == 323 || // K*+-
               pdg == 130 || // KL0
               pdg == 310 || // KS0
               pdg == 311 || // K0
               pdg == 313 || // K*0
               abs(pdg) == 221 || // eta
               abs(pdg) == 331 // eta' (958)
               ) {
        genie_n_mesons++;
      }
      if(pdg==211){
        nPip++;
      }

    }

    if (nPip != 1 || genie_n_mesons!= 1)
      return false;
    if (genie_n_photons != 0 )
      return false;

    int TargetPDG = kTruth_Target(sr);
    int TargetA = ((TargetPDG % 10000)) / 10;
    if(TargetA==1) return false;

    return true;


  }

  const TruthVar kTruth_IsSPP([](const caf::SRTrueInteractionProxy *nu) -> int {
    bool isspp = IsSPP(nu);
    if(isspp) return 1;
    else return 0;
  });

  const Var kNuMITrueIsSPP([](const caf::SRSliceProxy* slc) -> int {
    return kTruth_IsSPP(&slc->truth);
  });


  double GetSPPQ2Reweight(double Q2_GeV2){

    double X = Q2_GeV2;
    if(Q2_GeV2>=3.0) X = 2.5;

    double this_rw = 1.;
    if( X < 0.025000) this_rw = 1.253255;
    else if( X >= 0.025000 && X < 0.050000) this_rw = 1.589738;
    else if( X >= 0.050000 && X < 0.100000) this_rw = 1.733869;
    else if( X >= 0.100000 && X < 0.200000) this_rw = 1.651728;
    else if( X >= 0.200000 && X < 0.300000) this_rw = 1.659705;
    else if( X >= 0.300000 && X < 0.400000) this_rw = 1.584229;
    else if( X >= 0.400000 && X < 0.500000) this_rw = 1.703793;
    else if( X >= 0.500000 && X < 0.700000) this_rw = 1.475510;
    else if( X >= 0.700000 && X < 1.000000) this_rw = 1.456727;
    else if( X >= 1.000000 && X < 1.300000) this_rw = 1.252215;
    else if( X >= 1.300000 && X < 2.000000) this_rw = 1.048199;
    else if( X >= 2.000000 && X < 3.000000) this_rw = 1.650489;
    else{
      this_rw = 1.650489;
    }

    return this_rw;

  }
  double GetSPPLowQ2Suppression(double Q2_GeV2, double sigma){

    static double const Q2_Max = 0.7;
    static double const Q2_t1 = 0;
    static double const Q2_t2 = 0.35;
    static double const R1 = 0.3;
    static double const R2 = 0.6;

    if ((Q2_GeV2 > Q2_Max) || (Q2_GeV2 < 0)) {
      return 1.;
    }

    double RQ2 = (R2 * ((Q2_GeV2 - Q2_t1) * (Q2_GeV2 - Q2_Max)) /
           ((Q2_t2 - Q2_t1) * (Q2_t2 - Q2_Max))) +
          (((Q2_GeV2 - Q2_t1) * (Q2_GeV2 - Q2_t2)) /
           ((Q2_Max - Q2_t1) * (Q2_Max - Q2_t2)));
    return 1. - sigma * ((1. - R1) * pow((1. - RQ2), 2));

  }

  double GetSPPTpiCHLinearFitReweight(double Tpi_GeV){

    // CH result

    static double const P0 = 1.319098;
    //static double const P0Err = 0.120327;
    static double const P1 = -2.743935;
    //static double const P1Err = 0.488502;

    double X = Tpi_GeV;
    if(X>=0.350) X = 0.350;

    double this_rw = P0 + P1 * X;
    if(this_rw<0) this_rw = 1.;

    return this_rw;

  }
  double GetSPPTpiFeLinearFitReweight(double Tpi_GeV){

    // Fe result

    static double const P0 = 1.293700;
    //static double const P0Err = 0.147533;
    static double const P1 = -2.675087;
    //static double const P1Err = 0.520483;

    double X = Tpi_GeV;
    if(X>=0.350) X = 0.350;

    double this_rw = P0 + P1 * X;
    if(this_rw<0) this_rw = 1.;

    return this_rw;

  }
  double GetSPPTpiPbLinearFitReweight(double Tpi_GeV){

    // Pb result

    static double const P0 = 0.63527749;
    //static double const P0Err = 0.163223;
    static double const P1 = +1.274898;
    //static double const P1Err = 1.367774;

    double X = Tpi_GeV;
    if(X>=0.350) X = 0.350;

    double this_rw = P0 + P1 * X;
    if(this_rw<0) this_rw = 1.;

    return this_rw;

  }
  double GetSPPTpiMINERvATemplateReweight(double Tpi_GeV){

    double X = Tpi_GeV*1000.; // GeV to MeV

    double this_rw = 1.;
    if( X < 10.000000) this_rw = 0.267183;
    else if( X >= 10.000000 && X < 15.000000) this_rw = 0.218322;
    else if( X >= 15.000000 && X < 20.000000) this_rw = 0.372796;
    else if( X >= 20.000000 && X < 25.000000) this_rw = 0.587210;
    else if( X >= 25.000000 && X < 30.000000) this_rw = 0.767524;
    else if( X >= 30.000000 && X < 36.000000) this_rw = 0.880305;
    else if( X >= 36.000000 && X < 42.000000) this_rw = 0.669767;
    else if( X >= 42.000000 && X < 48.000000) this_rw = 0.817111;
    else if( X >= 48.000000 && X < 54.000000) this_rw = 1.092730;
    else if( X >= 54.000000 && X < 60.000000) this_rw = 0.995627;
    else if( X >= 60.000000 && X < 66.000000) this_rw = 0.916708;
    else if( X >= 66.000000 && X < 72.000000) this_rw = 1.243540;
    else if( X >= 72.000000 && X < 78.000000) this_rw = 1.211460;
    else if( X >= 78.000000 && X < 84.000000) this_rw = 1.121870;
    else if( X >= 84.000000 && X < 90.000000) this_rw = 1.253250;
    else if( X >= 90.000000 && X < 96.000000) this_rw = 1.191510;
    else if( X >= 96.000000 && X < 102.000000) this_rw = 1.038230;
    else if( X >= 102.000000 && X < 110.000000) this_rw = 1.237920;
    else if( X >= 110.000000 && X < 125.000000) this_rw = 1.190560;
    else if( X >= 125.000000 && X < 140.000000) this_rw = 1.229080;
    else if( X >= 140.000000 && X < 155.000000) this_rw = 0.988201;
    else if( X >= 155.000000 && X < 175.000000) this_rw = 1.032940;
    else if( X >= 175.000000 && X < 200.000000) this_rw = 0.901374;
    else if( X >= 200.000000 && X < 225.000000) this_rw = 0.757748;
    else if( X >= 225.000000 && X < 250.000000) this_rw = 0.755932;
    else if( X >= 250.000000 && X < 275.000000) this_rw = 0.638574;
    else if( X >= 275.000000 && X < 300.000000) this_rw = 0.493987;
    else if( X >= 300.000000 && X < 325.000000) this_rw = 0.391947;
    else if( X >= 325.000000 && X < 350.000000) this_rw = 0.323265;
    else if( X >= 350.000000 && X < 400.000000) this_rw = 0.452765;
    else if( X >= 400.000000 && X < 500.000000) this_rw = 0.594541;
    else if( X >= 500.000000 && X < 700.000000) this_rw = 0.768459;
    else if( X >= 700.000000 && X < 1000.000000) this_rw = 0.658024;
    else this_rw = 0.873622;

    return this_rw;

  }

  double GetSPPTpiMINERvAFittedReweight(double Tpi_GeV){

    static double landau_Cutoff = 0.225;

    if(Tpi_GeV<landau_Cutoff){
      // Params for Function = norm * ROOT.TMath.Landau(value, mu, sigma)
      // norm, mpv, width
      static double LandauParams[3] = {6.70797696, 0.12235454, 0.05731087};
      return LandauParams[0] * TMath::Landau(Tpi_GeV, LandauParams[1], LandauParams[2]);
    }
    else{
      if( landau_Cutoff <= Tpi_GeV && Tpi_GeV < 0.250000 ) return 0.755932;
      else if( 0.250000 <= Tpi_GeV && Tpi_GeV < 0.275000 ) return 0.638574;
      else if( 0.275000 <= Tpi_GeV && Tpi_GeV < 0.300000 ) return 0.493987;
      else if( 0.300000 <= Tpi_GeV && Tpi_GeV < 0.325000 ) return 0.391947;
      else if( 0.325000 <= Tpi_GeV && Tpi_GeV < 0.350000 ) return 0.323265;
      else if( 0.350000 <= Tpi_GeV && Tpi_GeV < 0.400000 ) return 0.452765;
      else if( 0.400000 <= Tpi_GeV && Tpi_GeV < 0.500000 ) return 0.594541;
      else if( 0.500000 <= Tpi_GeV && Tpi_GeV < 0.700000 ) return 0.768459;
      else if( 0.700000 <= Tpi_GeV && Tpi_GeV < 1.000000 ) return 0.658024;
      else if( 1.000000 <= Tpi_GeV && Tpi_GeV < 2.000000 ) return 0.873622;
      else return 0.873622;
    }

  }


  NuMIXSecPiSyst::NuMIXSecPiSyst(const std::string& name, const std::string& latexName):
    ISyst(name, latexName)
  {

  }

  void NuMIXSecPiSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    this->Shift(sigma, &sr->truth, weight);
  }

  void NuMIXSecPiSyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

    if( !IsSPP(sr) ) return;

    double CVCorr = kTruth_NuMISPPCVCorrection(sr);

    // 1/CVCorr is the correction back to nominal = 1sigma 
    double oneSigRW = 1./CVCorr;
    // Size of the one-sigma uncertainty obtained by subtracting 1
    double oneSigUnc = oneSigRW-1.;

    double this_rw = 1. + sigma * oneSigUnc;

    weight *= this_rw;

  }

  NuMIXSecLowQ2Suppression::NuMIXSecLowQ2Suppression(const std::string& name, const std::string& latexName):
    ISyst(name, latexName)
  {

  }

  void NuMIXSecLowQ2Suppression::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    this->Shift(sigma, &sr->truth, weight);
  }

  void NuMIXSecLowQ2Suppression::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

    if( !IsSPP(sr) ) return;

    double CVCorr = kTruth_NuMISPPLowQ2Suppression(sr);
    double Suppression = 1. - CVCorr;

    double this_rw = 1. - sigma * Suppression;

    weight *= this_rw;

  }


  // CV correction
  const TruthVar kTruth_NuMISPPCVCorrection([](const caf::SRTrueInteractionProxy *nu) -> float {

    if( !IsSPP(nu) ) return 1.;

    // Q2
    double Q2 = kTruth_Q2(nu);
    double Q2RW = GetSPPQ2Reweight(Q2);

    // Tpi
    double Tpi = kTruth_ChargedPionKE(nu);
    double TpiRW = GetSPPTpiMINERvAFittedReweight(Tpi);

    // CV correction as the full RW
    double CVCorr = Q2RW * TpiRW;

    return CVCorr;

  });
  const Var kNuMISPPCVCorrection([](const caf::SRSliceProxy* slc) -> float {

    return kTruth_NuMISPPCVCorrection(&slc->truth);

  });

  const TruthVar kTruth_NuMISPPLowQ2Suppression([](const caf::SRTrueInteractionProxy *nu) -> float {

    if( !IsSPP(nu) ) return 1.;

    // Q2
    double Q2 = kTruth_Q2(nu);
    return GetSPPLowQ2Suppression(Q2, 1.);

  });
  const Var kNuMISPPLowQ2Suppression([](const caf::SRSliceProxy* slc) -> float {

    return kTruth_NuMISPPLowQ2Suppression(&slc->truth);

  });

  

  // Separate reweights for study

  const Var kNuMISPPQ2RW([](const caf::SRSliceProxy* slc) -> float {

    if( !IsSPP(&slc->truth) ) return 1.;

    double Q2 = kTruth_Q2(&slc->truth);

    double Q2RW = GetSPPQ2Reweight(Q2); // Use Q2 from CH

    return Q2RW;

  });

  const Var kNuMISPPTpiCHLinearFitReweight([](const caf::SRSliceProxy* slc) -> float {

    if( !IsSPP(&slc->truth) ) return 1.;

    double Tpi = kTruth_ChargedPionKE(&slc->truth);

    double TpiRW = GetSPPTpiCHLinearFitReweight(Tpi);

    return TpiRW;

  });

  const Var kNuMISPPTpiFeLinearFitReweight([](const caf::SRSliceProxy* slc) -> float {

    if( !IsSPP(&slc->truth) ) return 1.;

    double Tpi = kTruth_ChargedPionKE(&slc->truth);

    double TpiRW = GetSPPTpiFeLinearFitReweight(Tpi);

    return TpiRW;

  });

  const Var kNuMISPPTpiPbLinearFitReweight([](const caf::SRSliceProxy* slc) -> float {

    if( !IsSPP(&slc->truth) ) return 1.;

    double Tpi = kTruth_ChargedPionKE(&slc->truth);

    double TpiRW = GetSPPTpiPbLinearFitReweight(Tpi);

    return TpiRW;

  });

  const Var kNuMISPPTpiMINERvATemplateReweight([](const caf::SRSliceProxy* slc) -> float {

    if( !IsSPP(&slc->truth) ) return 1.;

    double Tpi = kTruth_ChargedPionKE(&slc->truth);

    double TpiRW = GetSPPTpiMINERvATemplateReweight(Tpi);

    return TpiRW;

  });
  const Var kNuMISPPTpiMINERvAFittedReweight([](const caf::SRSliceProxy* slc) -> float {

    if( !IsSPP(&slc->truth) ) return 1.;

    double Tpi = kTruth_ChargedPionKE(&slc->truth);

    double TpiRW = GetSPPTpiMINERvAFittedReweight(Tpi);

    return TpiRW;

  });

  //---------------------------------------------------------------------
  // Split-track reweighting

  NuMIXSecSplitTrackReweight::NuMIXSecSplitTrackReweight(){
    // TH1* fRWCathode[2]; // [cryo]

    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIXSecSplitTrackReweight: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort();
    }

    // cathode
    std::string fCathodeFilePath = std::string(sbndata) +
                   "anaData/NuMI/TrackSplitReweight_cathode.root";
    TFile fCathode(fCathodeFilePath.c_str());
    if (fCathode.IsZombie()) {
      std::cout << "NuMIXSecSplitTrackReweight: Failed to open " << fCathodeFilePath << std::endl;
      std::abort();
    }

    for (int i_cryo : {0, 1}) {
      std::string strCryo = i_cryo==0 ? "East" : "West";
      for (int isSplit : {0, 1}) {
        std::string strIsSplit = isSplit==1 ? "Split" : "Reco";

        std::string hname = "RW_"+strIsSplit+"_"+strCryo;

        TH1* h = (TH1*)fCathode.Get(hname.c_str());
        if (!h) {
          std::cout << "NuMIXSecSplitTrackReweight: failed to find " << hname << " in " << fCathode.GetName()
                    << std::endl;
          std::abort();
        }
        h = (TH1*)h->Clone(UniqueName().c_str());
        h->SetDirectory(0);

        fRWCathode[i_cryo][isSplit] = h;

      }
    }

    // z=0
    std::string fZZeroFilePath = std::string(sbndata) +
                   "anaData/NuMI/TrackSplitReweight_zzero.root";
    TFile fZZero(fZZeroFilePath.c_str());
    if (fZZero.IsZombie()) {
      std::cout << "NuMIXSecSplitTrackReweight: Failed to open " << fZZeroFilePath << std::endl;
      std::abort();
    }

    for (int i_cryo : {0, 1}) {
      std::string strCryo = i_cryo==0 ? "East" : "West";
      for (int isSplit : {0, 1}) {
        std::string strIsSplit = isSplit==1 ? "Split" : "Reco";

        std::string hname = "RW_"+strIsSplit+"_"+strCryo;

        TH1* h = (TH1*)fZZero.Get(hname.c_str());
        if (!h) {
          std::cout << "NuMIXSecSplitTrackReweight: failed to find " << hname << " in " << fZZero.GetName()
                    << std::endl;
          std::abort();
        }
        h = (TH1*)h->Clone(UniqueName().c_str());
        h->SetDirectory(0);

        fRWZZero[i_cryo][isSplit] = h;

      }
    }

  }
  NuMIXSecSplitTrackReweight::~NuMIXSecSplitTrackReweight(){

  }
  NuMIXSecSplitTrackReweight& NuMIXSecSplitTrackReweight::Instance()
  {
    static NuMIXSecSplitTrackReweight splitTrackRW;
    return splitTrackRW;
  }
  int NuMIXSecSplitTrackReweight::CathodeSplitType(const caf::Proxy<caf::SRTrack>& trk) const{

    // -1: reweight not applicable
    //  0: reco-ed
    //  1: split

    double start_x = trk.start.x;
    double end_x = trk.end.x;
    double cathode_x = start_x>0. ? 210.21500 : -210.21500;

    bool CathodeCrossing = (start_x-cathode_x) * (end_x-cathode_x) < 0;
    bool EndAtCathode = abs(end_x)>207. && abs(end_x)<212.;

    // If cathode crossing, it is not split but reco-ed
    bool IsReco = CathodeCrossing;
    // If not crossing and end point is close to the cathode, it's a split track
    bool IsSplit = !CathodeCrossing && EndAtCathode;

    if(!IsReco && !IsSplit){
      return -1;
    }
    if(IsReco) return 0;
    if(IsSplit) return 1;

    // should not happen
    std::cout << "[NuMIXSecSplitTrackReweight::CathodeSplitType] Wrong track type" << std::endl;
    abort();
    return -2;

  }
  double NuMIXSecSplitTrackReweight::GetCathodeRW(const caf::Proxy<caf::SRTrack>& trk) const{

    int splitType = CathodeSplitType(trk);

    if(splitType<0){
      return 1.;
    }

    double cthetax = trk.dir.x;

    // x<0: East cryo = 0, x>0: West cryo = 1
    int idx_cryo = trk.start.x<0. ? 0 : 1;
    // Reco:0, Split:1
    int idx_IsSplit = splitType==1? 1 : 0;

    int this_bin = fRWCathode[idx_cryo][idx_IsSplit]->FindBin(cthetax);
    double rw = fRWCathode[idx_cryo][idx_IsSplit]->GetBinContent(this_bin);

    return rw;

  }

  int NuMIXSecSplitTrackReweight::ZZeroSplitType(const caf::Proxy<caf::SRTrack>& trk) const{

    // -1: reweight not applicable
    //  0: reco-ed
    //  1: split

    double start_z = trk.start.z;
    double end_z = trk.end.z;

    bool ZZeroCrossing = (start_z) * (end_z) < 0;
    bool EndAtZZero = end_z>-4. && end_z<4.;

    // If z=0 crossing, it is not split but reco-ed
    bool IsReco = ZZeroCrossing;
    // If not crossing and end point is close to z=0, it's a split track
    bool IsSplit = !ZZeroCrossing && EndAtZZero;

    if(!IsReco && !IsSplit){
      return -1;
    }
    if(IsReco) return 0;
    if(IsSplit) return 1;

    // should not happen
    std::cout << "[NuMIXSecSplitTrackReweight::ZZeroSplitType] Wrong track type" << std::endl;
    abort();
    return -2;

  }
  double NuMIXSecSplitTrackReweight::GetZZeroRW(const caf::Proxy<caf::SRTrack>& trk) const{

    int splitType = ZZeroSplitType(trk);

    if(splitType<0){
      return 1.;
    }
    // Z=0 boundary is not simulated in the MC
    // We only reweight the reco-ed track in the MC to match the data
    if(splitType==1){
      return 1.;
    }

    double cthetaz = trk.dir.z;

    // x<0: East cryo = 0, x>0: West cryo = 1
    int idx_cryo = trk.start.x<0. ? 0 : 1;
    // Reco:0, Split:1
    int idx_IsSplit = splitType==1? 1 : 0;

    int this_bin = fRWZZero[idx_cryo][idx_IsSplit]->FindBin(cthetaz);
    double rw = fRWZZero[idx_cryo][idx_IsSplit]->GetBinContent(this_bin);

    return rw;

  }

  // CV correction
  const Var kNuMISplitTrackCVCorrection([](const caf::SRSliceProxy* slc) -> float {
    int MuonIdx = kNuMIMuonCandidateIdx(slc);
    if(MuonIdx<0) return 1.;
    auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;

    const NuMIXSecSplitTrackReweight& splitTrackRW  = NuMIXSecSplitTrackReweight::Instance();

    //int cathodeSplitType = splitTrackRW.CathodeSplitType(trk);
    double rwCathodeSplit = splitTrackRW.GetCathodeRW(trk);
    //int zzeroSplitType = splitTrackRW.ZZeroSplitType(trk);
    double rwZZeroSplit = splitTrackRW.GetZZeroRW(trk);

    return rwCathodeSplit * rwZZeroSplit;

/*
    // When Z-crossing
    if(zzeroSplitType==0){

      // Also cathode crossing
      if(cathodeSplitType==0){

      }
      // But split at the cathode 
      if(cathodeSplitType==1){

      }

    }
*/

  });
  // Cathode and Z separately
  const Var kNuMICathodeSplitTrackCVCorrection([](const caf::SRSliceProxy* slc) -> float {
    int MuonIdx = kNuMIMuonCandidateIdx(slc);
    if(MuonIdx<0) return 1.;
    auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;

    const NuMIXSecSplitTrackReweight& splitTrackRW  = NuMIXSecSplitTrackReweight::Instance();

    return splitTrackRW.GetCathodeRW(trk);

  });
  const Var kNuMIZZeroSplitTrackCVCorrection([](const caf::SRSliceProxy* slc) -> float {
    int MuonIdx = kNuMIMuonCandidateIdx(slc);
    if(MuonIdx<0) return 1.;
    auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;

    const NuMIXSecSplitTrackReweight& splitTrackRW  = NuMIXSecSplitTrackReweight::Instance();

    return splitTrackRW.GetZZeroRW(trk);
    
  });


/*
  NuMIXSecSplitTrackSyst::NuMIXSecSplitTrackSyst(const std::string& name, const std::string& latexName):
    ISyst(name, latexName)
  {

  }

  void NuMIXSecSplitTrackSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {

    int MuonIdx = kNuMIMuonCandidateIdx(sr);
    if(MuonIdx<0) return;

    auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;


  }

  void NuMIXSecSplitTrackSyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {
  }
*/

} // end namespace ana
