// Bruce Howard 2024 - based on Calorimetry Syst (and others...)

/*
  ** NOTES during dev:
  TODO (26 Feb 2024): Switch to use the NORMAL vector with respect to the plane in question, and angle w.r.t. this, so for cathode basically dot product with the x direction...

  ** DESCRIPTION:
  TO UTLIZE THIS, ONE SHOULD RUN CAFMAKER IN A MODE THAT SAVES ALL THE CALO 
  POINTS TO THE TRACK CALO. THESE ARE USED THROUGHOUT TO PERFORM MODIFICATIONS.

  TODO: make a second track from the remnants.

  Some notes on what we do or do not modify:

  0. If the track crosses the cathode or z=0 point in CV reco, then throw a random 
     number and compare against the p for splitting in (z,theta_cathode) for the 
     cathode or (x, theta_(z=0)) for the z=0 gap.
        a. If it fails the comparison, do nothing
        b. If it passes the comparison, we split the track. Pick the calo point 
           closest to the cathode on each plane (TODO: use cathode bending maps to 
           pick point closest to post-bending instead) and this will be the new
           end of the track.

  1. We modify some reconstructed track quantities in SRTrack directly:
        producer = UNCHANGED
        npts = Use average number of calo points between ind2 and collection
	// CORRECTION 19 March 2024: trying SUM of all planes ncalopts... See if that matches better in the calc test...
        len = For each plane is [rr(first) - rr(final)] calo point. Avg ind2, coll.
        costh = UNCHANGED
        phi = UNCHANGED
        dir = UNCHANGED
        dir_end = set this to NaN for now
        start = UNCHANGED
        end = Average ind2 and collection for x, y, z points from final calo point

        bestplane = Update to plane with most calo nhits (see below)

  2. For collections that are saved to the SRTrack:
        SRTrkChi2PID = Use same method as Calo Syst (it should only use last 26cm of calo points)
        SRTrackCalo = Save new set of calo points. Update nhit n points with dedx < 1000.
            Set ke and charge to NaN for now (TODO: also update these...)
            See https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/CAFMaker/FillReco.cxx#L770
        SRTrkMCS = Set values to NaN for now
        SRTrkRange = Update this based on the new lengths

        Leave https://github.com/SBNSoftware/sbnanaobj/blob/develop/sbnanaobj/StandardRecord/SRTrack.h#L54-L62
        either unchanged or set to NaN for now:
            SRTrackTruth = UNCHANGED [[NOTE: match quantities will be wrong... hopefully true particle still relevant]]
            SRCRTHitMatch, SRCRTTrackMatch = UNCHANGED [[FOR NOW]]
            SRTrackScatterClosestApproach, SRTrackStoppingChi2Fit, SRTrackDazzle = Values set to NaN [[FOR NOW]]

TODO:
  4. See if we can make a second track out of the remaining bits. Here we will just edit the following things:
        STATUS: got it to the point of saving the info needed to make the below, so now need to do that...

        slc.reco.npfp += 1

        NOTE: We're keeping the same pfp-level things for the split track, like the track score & etc... That's maybe
        not ideal but we can't really recompute (at least at the moment)... For now we'll assume that both the split track
        and the separated track have the same pfp score info as the initial joint track. Also:
          parent = initial track id
          is_parent_primary = false
          etc.

        SRTrkChi2PID = Use same method as Calo Syst (it should only use last 26cm of calo points)
        SRTrackCalo = Save new set of calo points. Update nhit n points with dedx < 1000.
            Set ke and charge to NaN for now (TODO: also update these...)
            See https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/CAFMaker/FillReco.cxx#L770
        SRTrkRange = Update this based on the new lengths

        producer = UNCHANGED from initial track
        npts = Use average number of calo points between ind2 and collection
        len = For each plane is [rr(first) - rr(final)] calo point. Avg ind2, coll.
        costh = Set to NaN for now...
        phi = Set to NaN for now...
        dir = Take point-to-point directions from ind2 and collection from first 5 points each, average.
        dir_end = UNCHANGED from initial track
        start = Average ind2 and collection for x, y, z points from final calo point
        end = UNCHANGED from initial track

        bestplane = Update to plane with most calo nhits (as above.)

  5. ALSO will want to deal with showers for how things will play into shower cut(s)...
    -- vtx to start point as conv gap
    -- 90% of length as shw length
    -- FOR THE SECOND PFP, setting the start point to be the start point of the second pfp's track (obvious bug fixed on 19 March 2024)
*/

#include "sbnana/CAFAna/Systs/TrackSplitSyst.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "cetlib/search_path.h"
#include <cmath>
#include <iostream>
#include "TMath.h"
#include <memory>

#include <array>
#include "TGraph.h"

#include "TVector3.h"

#include "sbnanaobj/StandardRecord/SRConstants.h"
#include "sbnanaobj/StandardRecord/SRTrackDazzle.h"
#include "sbnanaobj/StandardRecord/SRTrkMCS.h"
#include "sbnanaobj/StandardRecord/SRTrackScatterClosestApproach.h"
#include "sbnanaobj/StandardRecord/SRTrackStoppingChi2Fit.h"
#include "sbnanaobj/StandardRecord/SRCaloPoint.h"
#include "sbnanaobj/StandardRecord/SRTrueCaloPoint.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"
#include "sbnanaobj/StandardRecord/SRPFP.h"

namespace ana {

  TrackSplitSyst::TrackSplitSyst(const std::string& name, const std::string& latexName, const bool debugPrint, const std::string& fileName):
    ISyst(name, latexName), fDebug(debugPrint)
  {
    std::string useFileName;

    // Load in the histogram, as in the NuMI Flux Uncertainty
    if ( fileName.size()==0 ) {
      const char* sbndata = std::getenv("SBNDATA_DIR");
      if (!sbndata) {
        std::cout << "NuMIPpfxFluxWeight: $SBNDATA_DIR environment variable not set. Please setup "
                     "the sbndata product."
                  << std::endl;
        std::abort();
      }

      useFileName = std::string(sbndata) + "/recoData/output_histo_split_tracks_data.root";
    }
    else useFileName = fileName;

    TFile* file_tracksplitdata = TFile::Open(useFileName.c_str());

    splitProb[0][0] = (TH2D*)file_tracksplitdata->Get("z_cathodeangle_east");
    splitProb[0][1] = (TH2D*)file_tracksplitdata->Get("x_zgapangle_east");
    splitProb[1][0] = (TH2D*)file_tracksplitdata->Get("z_cathodeangle_west");
    splitProb[1][1] = (TH2D*)file_tracksplitdata->Get("x_zgapangle_west");

    // TEMPORARY VERSION WHERE WE ASSUME A CONSTANT 15% SPLITTING AT CATHODE EVERYWHERE
    // For each cryo:
    //   20 bins along z for cathode (z=0 uses 10 bins across x)
    //   5 bins of angle
    /*
    splitProb[0][0] = new TH2D("TrackSplitting_CryoEast_Cathode", ";Z [cm];cosThCathode", 20, -900, 900, 5, -1, 1);
    splitProb[0][1] = new TH2D("TrackSplitting_CryoEast_MiddleZ", ";X [cm];cosThMiddleZ", 10, -360, -60, 5, -1, 1);
    splitProb[1][0] = new TH2D("TrackSplitting_CryoWest_Cathode", ";Z [cm];cosThCathode", 20, -900, 900, 5, -1, 1);
    splitProb[1][1] = new TH2D("TrackSplitting_CryoWest_MiddleZ", ";X [cm];cosThMiddleZ", 10, 60, 360, 5, -1, 1);
    for ( unsigned int idxX = 1; idxX <= 20; ++idxX ) {
      for ( unsigned int idxY = 1; idxY <= 10; ++idxY ) {
        splitProb[0][0]->SetBinContent(idxX, idxY, 0.15);
        splitProb[1][0]->SetBinContent(idxX, idxY, 0.15);
        if ( idxX <= 10 ) {
          splitProb[0][1]->SetBinContent(idxX, idxY, 0.15);
          splitProb[1][1]->SetBinContent(idxX, idxY, 0.15);
        }
      }
    }
    */

    // Get the dEdx templates as in Calo Syst
    cet::search_path sp("FW_SEARCH_PATH");

    std::string kChi2TemplateFileName = "dEdxrestemplates.root";
    std::string kChi2TemplateFullFilePath;
    sp.find_file(kChi2TemplateFileName, kChi2TemplateFullFilePath);

    TFile *file_Chi2Template = TFile::Open(kChi2TemplateFullFilePath.c_str());
    dedx_range_pro = (TProfile*)file_Chi2Template->Get("dedx_range_pro");
    dedx_range_ka  = (TProfile*)file_Chi2Template->Get("dedx_range_ka");
    dedx_range_pi  = (TProfile*)file_Chi2Template->Get("dedx_range_pi");
    dedx_range_mu  = (TProfile*)file_Chi2Template->Get("dedx_range_mu");

    // Set up the random engine
    tRand = new TRandom3(0);

    // Set up the spline used for momentum calculation, as in
    // https://github.com/LArSoft/larreco/blob/LARSOFT_SUITE_v09_72_00/larreco/RecoAlg/TrackMomentumCalculator.cxx
    std::array<float, 29> Range_grampercm{
      {9.833E-1/1.396, 1.786E0/1.396, 3.321E0/1.396, 6.598E0/1.396, 1.058E1/1.396, 3.084E1/1.396, 4.250E1/1.396, 6.732E1/1.396, 1.063E2/1.396, 1.725E2/1.396,
       2.385E2/1.396,  4.934E2/1.396, 6.163E2/1.396, 8.552E2/1.396, 1.202E3/1.396, 1.758E3/1.396, 2.297E3/1.396, 4.359E3/1.396, 5.354E3/1.396, 7.298E3/1.396,
       1.013E4/1.396,  1.469E4/1.396, 1.910E4/1.396, 3.558E4/1.396, 4.326E4/1.396, 5.768E4/1.396, 7.734E4/1.396, 1.060E5/1.396, 1.307E5/1.396}};

    constexpr std::array<float, 29> KE_MeV{
      {10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
      400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
      20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000}};
    TGraph const KEvsR{29, Range_grampercm.data(), KE_MeV.data()};
    KEvsR_spline3 = TSpline3("KEvsRS", &KEvsR);
  }

  TrkChi2Results TrackSplitSyst::CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo) const {
    // copied&modified from https://github.com/LArSoft/larana/blob/develop/larana/ParticleIdentification/Chi2PIDAlg.cxx#L60

    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double PIDA = 0; //by Bruce Baller
    std::vector<double> vpida;

    int used_trkres = 0;
    for(unsigned i = 0; i < calo.points.size(); ++i) { //hits
      const auto& pt = calo.points[i];
      double hit_dedx = pt.dedx;
      double hit_rr = pt.rr;

      //ignore the first and the last point
      if (i == 0 || i == calo.points.size() - 1) continue;
      if (hit_rr < 30) {
        PIDA += hit_dedx * pow(hit_rr, 0.42);
        vpida.push_back(hit_dedx * pow(hit_rr, 0.42));
        used_trkres++;
      }
      if (hit_dedx > 1000) continue; //protect against large pulse height
      int bin = dedx_range_pro->FindBin(hit_rr);
      if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
        double bincpro = dedx_range_pro->GetBinContent(bin);
        if (bincpro < 1e-6) { //for 0 bin content, using neighboring bins
          bincpro =
            (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
        }
        double bincka = dedx_range_ka->GetBinContent(bin);
        if (bincka < 1e-6) {
          bincka =
            (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;
        }
        double bincpi = dedx_range_pi->GetBinContent(bin);
        if (bincpi < 1e-6) {
          bincpi =
            (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
        }
        double bincmu = dedx_range_mu->GetBinContent(bin);
        if (bincmu < 1e-6) {
          bincmu =
            (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
        }
        double binepro = dedx_range_pro->GetBinError(bin);
        if (binepro < 1e-6) {
          binepro =
            (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;
        }
        double bineka = dedx_range_ka->GetBinError(bin);
        if (bineka < 1e-6) {
          bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;
        }
        double binepi = dedx_range_pi->GetBinError(bin);
        if (binepi < 1e-6) {
          binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;
        }
        double binemu = dedx_range_mu->GetBinError(bin);
        if (binemu < 1e-6) {
          binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;
        }
        //double errke = 0.05*hit_dedx;   //5% KE resolution
        double errdedx = 0.04231 + 0.0001783 * hit_dedx * hit_dedx; //resolution on dE/dx
        errdedx *= hit_dedx;
        chi2pro += pow((hit_dedx - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
        chi2ka += pow((hit_dedx - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
        chi2pi += pow((hit_dedx - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
        chi2mu += pow((hit_dedx - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
        //std::cout<<i<<" "<<hit_dedx<<" "<<hit_rr<<" "<<bincpro<<std::endl;
        ++npt;
      }
    }

    // Making output
    TrkChi2Results output;
    if (npt) {
      chi2pro /= npt;
      chi2ka /= npt;
      chi2pi /= npt;
      chi2mu /= npt;
    }
    static bool fUseMedian = true; // I think default fcl in ICARUS is true
    if (used_trkres > 0) {
      if (fUseMedian) {
        PIDA = TMath::Median(vpida.size(), &vpida[0]);
      }
      else { // use mean
        PIDA /= used_trkres;
      }
    }

    output.chi2_kaon = chi2ka;
    output.chi2_muon = chi2mu;
    output.chi2_pion = chi2pi;
    output.chi2_proton = chi2pro;
    output.pida = PIDA;
    output.pid_ndof = npt;

    return output;
  } // Calculate Chi2

  TrkChi2Results TrackSplitSyst::CalculateChi2(const caf::SRTrackCalo& calo) const {
    // copied&modified from https://github.com/LArSoft/larana/blob/develop/larana/ParticleIdentification/Chi2PIDAlg.cxx#L60

    int npt = 0;
    double chi2pro = 0;
    double chi2ka = 0;
    double chi2pi = 0;
    double chi2mu = 0;
    double PIDA = 0; //by Bruce Baller
    std::vector<double> vpida;

    int used_trkres = 0;
    for(unsigned i = 0; i < calo.points.size(); ++i) { //hits
      const auto& pt = calo.points[i];
      double hit_dedx = pt.dedx;
      double hit_rr = pt.rr;

      //ignore the first and the last point
      if (i == 0 || i == calo.points.size() - 1) continue;
      if (hit_rr < 30) {
        PIDA += hit_dedx * pow(hit_rr, 0.42);
        vpida.push_back(hit_dedx * pow(hit_rr, 0.42));
        used_trkres++;
      }
      if (hit_dedx > 1000) continue; //protect against large pulse height
      int bin = dedx_range_pro->FindBin(hit_rr);
      if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
        double bincpro = dedx_range_pro->GetBinContent(bin);
        if (bincpro < 1e-6) { //for 0 bin content, using neighboring bins
          bincpro =
            (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
        }
        double bincka = dedx_range_ka->GetBinContent(bin);
        if (bincka < 1e-6) {
          bincka =
            (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;
        }
        double bincpi = dedx_range_pi->GetBinContent(bin);
        if (bincpi < 1e-6) {
          bincpi =
            (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
        }
        double bincmu = dedx_range_mu->GetBinContent(bin);
        if (bincmu < 1e-6) {
          bincmu =
            (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
        }
        double binepro = dedx_range_pro->GetBinError(bin);
        if (binepro < 1e-6) {
          binepro =
            (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;
        }
        double bineka = dedx_range_ka->GetBinError(bin);
        if (bineka < 1e-6) {
          bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;
        }
        double binepi = dedx_range_pi->GetBinError(bin);
        if (binepi < 1e-6) {
          binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;
        }
        double binemu = dedx_range_mu->GetBinError(bin);
        if (binemu < 1e-6) {
          binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;
        }
        //double errke = 0.05*hit_dedx;   //5% KE resolution
        double errdedx = 0.04231 + 0.0001783 * hit_dedx * hit_dedx; //resolution on dE/dx
        errdedx *= hit_dedx;
        chi2pro += pow((hit_dedx - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
        chi2ka += pow((hit_dedx - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
        chi2pi += pow((hit_dedx - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
        chi2mu += pow((hit_dedx - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
        //std::cout<<i<<" "<<hit_dedx<<" "<<hit_rr<<" "<<bincpro<<std::endl;
        ++npt;
      }
    }

    // Making output
    TrkChi2Results output;
    if (npt) {
      chi2pro /= npt;
      chi2ka /= npt;
      chi2pi /= npt;
      chi2mu /= npt;
    }
    static bool fUseMedian = true; // I think default fcl in ICARUS is true
    if (used_trkres > 0) {
      if (fUseMedian) {
        PIDA = TMath::Median(vpida.size(), &vpida[0]);
      }
      else { // use mean
        PIDA /= used_trkres;
      }
    }

    output.chi2_kaon = chi2ka;
    output.chi2_muon = chi2mu;
    output.chi2_pion = chi2pi;
    output.chi2_proton = chi2pro;
    output.pida = PIDA;
    output.pid_ndof = npt;

    return output;
  } // Calculate Chi2

  TrkMomentumResults TrackSplitSyst::CalculateMomenta(const float length) const {
    // largely copying and modifying from
    // https://github.com/LArSoft/larreco/blob/develop/larreco/RecoAlg/TrackMomentumCalculator.cxx
    // from LARSOFT_SUITE_v09_72_00. The track momentum calculator is used by:
    // https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/LArRecoProducer/RangePAllPID_module.cc
    // and we copy its implementation of pion momentum (the larsoft module covers muon and proton...)

    // Re-implementing pieces to store the values, but the calculation bits should be the same.

    if (length < 0 || std::isnan(length)) {
      std::cout << "Invalid track range " << length << " return -1" << std::endl;
      TrkMomentumResults output;
      output.p_proton = -1.;
      output.p_pion = -1.;
      output.p_muon = -1.;

      return output;
    }

    constexpr double Muon_M = 105.7, Proton_M = 938.272;
    constexpr double Pion_M = 139.57; // https://pdg.lbl.gov/2018/listings/rpp2018-list-pi-plus-minus.pdf

    TrkMomentumResults output;

    // Muon
    double KE_Muon = KEvsR_spline3.Eval(length);
    double p_muon = -999.;
    if ( KE_Muon >= 0 ) p_muon = std::sqrt((KE_Muon * KE_Muon) + (2 * Muon_M * KE_Muon));
    output.p_muon = p_muon / 1000;

    // Pion, as implemented in sbncode RangeP code
    double KE_Pion_FromMuCalc = KEvsR_spline3.Eval(length*Muon_M/Pion_M);
    double p_pion = -999.;
    if ( KE_Pion_FromMuCalc >= 0 )
      p_pion = std::sqrt((KE_Pion_FromMuCalc * KE_Pion_FromMuCalc) + (2 * Muon_M * KE_Pion_FromMuCalc)) * (Pion_M/Muon_M);
    output.p_pion = p_pion / 1000;

    // Proton
    double KE_Proton = -999.;
    if (length > 0 && length <= 80)
      KE_Proton = 29.9317 * std::pow(length, 0.586304);
    else if (length > 80 && length <= 3.022E3)
      KE_Proton = 149.904 + (3.34146 * length) + (-0.00318856 * length * length) +
                  (4.34587E-6 * length * length * length) +
                  (-3.18146E-9 * length * length * length * length) +
                  (1.17854E-12 * length * length * length * length * length) +
                  (-1.71763E-16 * length * length * length * length * length * length);
    double p_proton = -999.;
    if ( KE_Proton >= 0 ) p_proton = std::sqrt((KE_Proton * KE_Proton) + (2 * Proton_M * KE_Proton));
    output.p_proton = p_proton / 1000;

    return output;
  }

  caf::SRCaloPoint TrackSplitSyst::FillCaloPointFrom( const caf::Proxy<caf::SRCaloPoint>& inCaloPt ) const {
    caf::SRCaloPoint outCaloPt;
    outCaloPt.rr = inCaloPt.rr;
    outCaloPt.dqdx = inCaloPt.dqdx;
    outCaloPt.dedx = inCaloPt.dedx;
    outCaloPt.pitch = inCaloPt.pitch;
    outCaloPt.t = inCaloPt.t;
    outCaloPt.x = inCaloPt.x;
    outCaloPt.y = inCaloPt.y;
    outCaloPt.z = inCaloPt.z;
    outCaloPt.integral = inCaloPt.integral;
    outCaloPt.sumadc = inCaloPt.sumadc;
    outCaloPt.width = inCaloPt.width;
    outCaloPt.mult = inCaloPt.mult;
    outCaloPt.wire = inCaloPt.wire;
    outCaloPt.tpc = inCaloPt.tpc;
    outCaloPt.start = inCaloPt.start;
    outCaloPt.end = inCaloPt.end;
    outCaloPt.channel = inCaloPt.channel;

    caf::SRTrueCaloPoint outTruth;
    outTruth.h_nelec = inCaloPt.truth.h_nelec;
    outTruth.h_e = inCaloPt.truth.h_nelec;
    outTruth.p_nelec = inCaloPt.truth.h_nelec;
    outTruth.h_e = inCaloPt.truth.h_nelec;
    outTruth.x = inCaloPt.truth.h_nelec;
    outTruth.y = inCaloPt.truth.h_nelec;
    outTruth.z = inCaloPt.truth.h_nelec;
    outTruth.rr = inCaloPt.truth.h_nelec;
    outTruth.pitch = inCaloPt.truth.h_nelec;

    outCaloPt.truth = outTruth;

    return outCaloPt;
  }

  void TrackSplitSyst::FillPtrPFP( caf::SRPFP& ret, const caf::Proxy<caf::SRPFP>& inPfp ) const
  {
    ret.id = inPfp.id;
    ret.ndaughters = inPfp.ndaughters;
    // Set ret.daughters = inPfp.daughters;
    for ( auto const& daughter : inPfp.daughters ) ret.daughters.push_back( int(daughter) );

    ret.parent = inPfp.parent;
    ret.parent_is_primary = inPfp.parent_is_primary;

    ret.trackScore = inPfp.trackScore;
    // Set ret.pfochar = inPfp.pfochar;
    ret.pfochar.chgendfrac = inPfp.pfochar.chgendfrac;
    ret.pfochar.chgfracspread = inPfp.pfochar.chgfracspread;
    ret.pfochar.linfitdiff = inPfp.pfochar.linfitdiff;
    ret.pfochar.linfitlen = inPfp.pfochar.linfitlen;
    ret.pfochar.linfitgaplen = inPfp.pfochar.linfitgaplen;
    ret.pfochar.linfitrms = inPfp.pfochar.linfitrms;
    ret.pfochar.openanglediff = inPfp.pfochar.openanglediff;
    ret.pfochar.pca2ratio = inPfp.pfochar.pca2ratio;
    ret.pfochar.pca3ratio = inPfp.pfochar.pca3ratio;
    ret.pfochar.vtxdist = inPfp.pfochar.vtxdist;

    ret.slcID = inPfp.slcID;
    ret.t0 = inPfp.t0;

    // Set ret.trk = inPfp.trk;
    ret.trk.producer = inPfp.trk.producer;
    ret.trk.npts = inPfp.trk.npts;
    ret.trk.len = inPfp.trk.len;
    ret.trk.costh = inPfp.trk.costh;
    ret.trk.phi = inPfp.trk.phi;
    ret.trk.dir.x = inPfp.trk.dir.x;
    ret.trk.dir.y = inPfp.trk.dir.y;
    ret.trk.dir.z = inPfp.trk.dir.z;
    ret.trk.dir_end.x = inPfp.trk.dir_end.x;
    ret.trk.dir_end.y = inPfp.trk.dir_end.y;
    ret.trk.dir_end.z = inPfp.trk.dir_end.z;
    ret.trk.start.x = inPfp.trk.start.x;
    ret.trk.start.y = inPfp.trk.start.y;
    ret.trk.start.z = inPfp.trk.start.z;
    ret.trk.end.x = inPfp.trk.end.x;
    ret.trk.end.y = inPfp.trk.end.y;
    ret.trk.end.z = inPfp.trk.end.z;
    ret.trk.bestplane = inPfp.trk.bestplane;
    // Trk Chi2PID
    ret.trk.chi2pid[0].pdg = inPfp.trk.chi2pid[0].pdg;
    ret.trk.chi2pid[0].pid_ndof = inPfp.trk.chi2pid[0].pid_ndof;
    ret.trk.chi2pid[0].chi2_muon = inPfp.trk.chi2pid[0].chi2_muon;
    ret.trk.chi2pid[0].chi2_pion = inPfp.trk.chi2pid[0].chi2_pion;
    ret.trk.chi2pid[0].chi2_kaon = inPfp.trk.chi2pid[0].chi2_kaon;
    ret.trk.chi2pid[0].chi2_proton = inPfp.trk.chi2pid[0].chi2_proton;
    ret.trk.chi2pid[0].pida = inPfp.trk.chi2pid[0].pida;
    ret.trk.chi2pid[1].pdg = inPfp.trk.chi2pid[1].pdg;
    ret.trk.chi2pid[1].pid_ndof = inPfp.trk.chi2pid[1].pid_ndof;
    ret.trk.chi2pid[1].chi2_muon = inPfp.trk.chi2pid[1].chi2_muon;
    ret.trk.chi2pid[1].chi2_pion = inPfp.trk.chi2pid[1].chi2_pion;
    ret.trk.chi2pid[1].chi2_kaon = inPfp.trk.chi2pid[1].chi2_kaon;
    ret.trk.chi2pid[1].chi2_proton = inPfp.trk.chi2pid[1].chi2_proton;
    ret.trk.chi2pid[1].pida = inPfp.trk.chi2pid[1].pida;
    ret.trk.chi2pid[2].pdg = inPfp.trk.chi2pid[2].pdg;
    ret.trk.chi2pid[2].pid_ndof = inPfp.trk.chi2pid[2].pid_ndof;
    ret.trk.chi2pid[2].chi2_muon = inPfp.trk.chi2pid[2].chi2_muon;
    ret.trk.chi2pid[2].chi2_pion = inPfp.trk.chi2pid[2].chi2_pion;
    ret.trk.chi2pid[2].chi2_kaon = inPfp.trk.chi2pid[2].chi2_kaon;
    ret.trk.chi2pid[2].chi2_proton = inPfp.trk.chi2pid[2].chi2_proton;
    ret.trk.chi2pid[2].pida = inPfp.trk.chi2pid[2].pida;
    // Trk Calo
    ret.trk.calo[0].nhit = inPfp.trk.calo[0].nhit;
    ret.trk.calo[0].ke = inPfp.trk.calo[0].ke;
    ret.trk.calo[0].charge = inPfp.trk.calo[0].charge;
    for ( auto const& caloPt : inPfp.trk.calo[0].points ) ret.trk.calo[0].points.push_back( FillCaloPointFrom(caloPt) );
    ret.trk.calo[1].nhit = inPfp.trk.calo[1].nhit;
    ret.trk.calo[1].ke = inPfp.trk.calo[1].ke;
    ret.trk.calo[1].charge = inPfp.trk.calo[1].charge;
    for ( auto const& caloPt : inPfp.trk.calo[1].points ) ret.trk.calo[1].points.push_back( FillCaloPointFrom(caloPt) );
    ret.trk.calo[2].nhit = inPfp.trk.calo[2].nhit;
    ret.trk.calo[2].ke = inPfp.trk.calo[2].ke;
    ret.trk.calo[2].charge = inPfp.trk.calo[2].charge;
    for ( auto const& caloPt : inPfp.trk.calo[2].points ) ret.trk.calo[2].points.push_back( FillCaloPointFrom(caloPt) );
    // MCS P
    ret.trk.mcsP.fwdP_muon = inPfp.trk.mcsP.fwdP_muon;
    ret.trk.mcsP.fwdP_pion = inPfp.trk.mcsP.fwdP_pion;
    ret.trk.mcsP.fwdP_kaon = inPfp.trk.mcsP.fwdP_kaon;
    ret.trk.mcsP.fwdP_proton = inPfp.trk.mcsP.fwdP_proton;
    ret.trk.mcsP.fwdP_err_muon = inPfp.trk.mcsP.fwdP_err_muon;
    ret.trk.mcsP.fwdP_err_pion = inPfp.trk.mcsP.fwdP_err_pion;
    ret.trk.mcsP.fwdP_err_kaon = inPfp.trk.mcsP.fwdP_err_kaon;
    ret.trk.mcsP.fwdP_err_proton = inPfp.trk.mcsP.fwdP_err_proton;
    ret.trk.mcsP.bwdP_muon = inPfp.trk.mcsP.bwdP_muon;
    ret.trk.mcsP.bwdP_pion = inPfp.trk.mcsP.bwdP_pion;
    ret.trk.mcsP.bwdP_kaon = inPfp.trk.mcsP.bwdP_kaon;
    ret.trk.mcsP.bwdP_proton = inPfp.trk.mcsP.bwdP_proton;
    ret.trk.mcsP.bwdP_err_muon = inPfp.trk.mcsP.bwdP_err_muon;
    ret.trk.mcsP.bwdP_err_pion = inPfp.trk.mcsP.bwdP_err_pion;
    ret.trk.mcsP.bwdP_err_kaon = inPfp.trk.mcsP.bwdP_err_kaon;
    ret.trk.mcsP.bwdP_err_proton = inPfp.trk.mcsP.bwdP_err_proton;
    ret.trk.mcsP.is_bwd_muon = inPfp.trk.mcsP.is_bwd_muon;
    ret.trk.mcsP.is_bwd_pion = inPfp.trk.mcsP.is_bwd_pion;
    ret.trk.mcsP.is_bwd_kaon = inPfp.trk.mcsP.is_bwd_kaon;
    ret.trk.mcsP.is_bwd_proton = inPfp.trk.mcsP.is_bwd_proton;
    // Range P
    ret.trk.rangeP.p_muon = inPfp.trk.rangeP.p_muon;
    ret.trk.rangeP.p_pion = inPfp.trk.rangeP.p_pion;
    ret.trk.rangeP.p_proton = inPfp.trk.rangeP.p_proton;
    // Trk Truth
    ret.trk.truth.visEintrk = inPfp.trk.truth.visEintrk;
    ret.trk.truth.eff = inPfp.trk.truth.eff;
    ret.trk.truth.eff_cryo = inPfp.trk.truth.eff_cryo;
    ret.trk.truth.pur = inPfp.trk.truth.pur;
    ret.trk.truth.nmatches = inPfp.trk.truth.nmatches;
    // --> Leave the "matches" empty...
    ret.trk.truth.bestmatch.G4ID = inPfp.trk.truth.bestmatch.G4ID;
    ret.trk.truth.bestmatch.energy = inPfp.trk.truth.bestmatch.energy;
    ret.trk.truth.bestmatch.hit_completeness = inPfp.trk.truth.bestmatch.hit_completeness;
    ret.trk.truth.bestmatch.hit_purity = inPfp.trk.truth.bestmatch.hit_purity;
    ret.trk.truth.bestmatch.energy_completeness = inPfp.trk.truth.bestmatch.energy_completeness;
    ret.trk.truth.bestmatch.energy_purity = inPfp.trk.truth.bestmatch.energy_purity;
    ret.trk.truth.p.genE = inPfp.trk.truth.p.genE;
    // --> Leave the plane info NaN'ed/zero'ed out...
    ret.trk.truth.p.startE = inPfp.trk.truth.p.startE;
    ret.trk.truth.p.endE = inPfp.trk.truth.p.endE;
    ret.trk.truth.p.genT = inPfp.trk.truth.p.genT;
    ret.trk.truth.p.startT = inPfp.trk.truth.p.startT;
    ret.trk.truth.p.endT = inPfp.trk.truth.p.endT;
    ret.trk.truth.p.length = inPfp.trk.truth.p.length;
    ret.trk.truth.p.genp.x = inPfp.trk.truth.p.genp.x;
    ret.trk.truth.p.genp.y = inPfp.trk.truth.p.genp.y;
    ret.trk.truth.p.genp.z = inPfp.trk.truth.p.genp.z;
    ret.trk.truth.p.startp.x = inPfp.trk.truth.p.startp.x;
    ret.trk.truth.p.startp.y = inPfp.trk.truth.p.startp.y;
    ret.trk.truth.p.startp.z = inPfp.trk.truth.p.startp.z;
    ret.trk.truth.p.endp.x = inPfp.trk.truth.p.endp.x;
    ret.trk.truth.p.endp.y = inPfp.trk.truth.p.endp.y;
    ret.trk.truth.p.endp.z = inPfp.trk.truth.p.endp.z;
    ret.trk.truth.p.gen.x = inPfp.trk.truth.p.gen.x;
    ret.trk.truth.p.gen.y = inPfp.trk.truth.p.gen.y;
    ret.trk.truth.p.gen.z = inPfp.trk.truth.p.gen.z;
    ret.trk.truth.p.start.x = inPfp.trk.truth.p.start.x;
    ret.trk.truth.p.start.y = inPfp.trk.truth.p.start.y;
    ret.trk.truth.p.start.z = inPfp.trk.truth.p.start.z;
    ret.trk.truth.p.end.x = inPfp.trk.truth.p.end.x;
    ret.trk.truth.p.end.y = inPfp.trk.truth.p.end.y;
    ret.trk.truth.p.end.z = inPfp.trk.truth.p.end.z;
    ret.trk.truth.p.wallin = inPfp.trk.truth.p.wallin;
    ret.trk.truth.p.wallout = inPfp.trk.truth.p.wallout;
    ret.trk.truth.p.cont_tpc = inPfp.trk.truth.p.cont_tpc;
    ret.trk.truth.p.crosses_tpc = inPfp.trk.truth.p.crosses_tpc;
    ret.trk.truth.p.contained = inPfp.trk.truth.p.contained;
    ret.trk.truth.p.pdg = inPfp.trk.truth.p.pdg;
    ret.trk.truth.p.G4ID = inPfp.trk.truth.p.G4ID;
    ret.trk.truth.p.interaction_id = inPfp.trk.truth.p.interaction_id;
    ret.trk.truth.p.cryostat = inPfp.trk.truth.p.cryostat;
    for ( auto const& daughter : inPfp.trk.truth.p.daughters ) ret.trk.truth.p.daughters.push_back( (unsigned)daughter );
    ret.trk.truth.p.parent = inPfp.trk.truth.p.parent;
    ret.trk.truth.p.generator = inPfp.trk.truth.p.generator;
    ret.trk.truth.p.start_process = inPfp.trk.truth.p.start_process;
    ret.trk.truth.p.end_process = inPfp.trk.truth.p.end_process;
    ret.trk.truth.p.gstatus = inPfp.trk.truth.p.gstatus;
    // --> For now leave other stuff NaN'ed/zero'ed out...

    // Set ret.shw = inPfp.shw;
    ret.shw.bestplane = inPfp.shw.bestplane;
    ret.shw.bestplane_dEdx = inPfp.shw.bestplane_dEdx;
    ret.shw.bestplane_energy = inPfp.shw.bestplane_energy;
    ret.shw.conversion_gap = inPfp.shw.conversion_gap;
    ret.shw.density = inPfp.shw.density;
    ret.shw.len = inPfp.shw.len;
    ret.shw.open_angle = inPfp.shw.open_angle;
    ret.shw.dir.x = inPfp.shw.dir.x;
    ret.shw.dir.y = inPfp.shw.dir.y;
    ret.shw.dir.z = inPfp.shw.dir.z;
    ret.shw.start.x = inPfp.shw.start.x;
    ret.shw.start.y = inPfp.shw.start.y;
    ret.shw.start.z = inPfp.shw.start.z;
    ret.shw.end.x = inPfp.shw.end.x;
    ret.shw.end.y = inPfp.shw.end.y;
    ret.shw.end.z = inPfp.shw.end.z;
    ret.shw.cosmicDist = inPfp.shw.cosmicDist;
    ret.shw.producer = inPfp.shw.producer;
    // Shw Truth
    ret.shw.truth.visEintrk = inPfp.shw.truth.visEintrk;
    ret.shw.truth.eff = inPfp.shw.truth.eff;
    ret.shw.truth.eff_cryo = inPfp.shw.truth.eff_cryo;
    ret.shw.truth.pur = inPfp.shw.truth.pur;
    ret.shw.truth.nmatches = inPfp.shw.truth.nmatches;
    // --> Leave the "matches" empty...
    ret.shw.truth.bestmatch.G4ID = inPfp.shw.truth.bestmatch.G4ID;
    ret.shw.truth.bestmatch.energy = inPfp.shw.truth.bestmatch.energy;
    ret.shw.truth.bestmatch.hit_completeness = inPfp.shw.truth.bestmatch.hit_completeness;
    ret.shw.truth.bestmatch.hit_purity = inPfp.shw.truth.bestmatch.hit_purity;
    ret.shw.truth.bestmatch.energy_completeness = inPfp.shw.truth.bestmatch.energy_completeness;
    ret.shw.truth.bestmatch.energy_purity = inPfp.shw.truth.bestmatch.energy_purity;
    ret.shw.truth.p.genE = inPfp.shw.truth.p.genE;
    // --> Leave the plane info NaN'ed/zero'ed out...
    ret.shw.truth.p.startE = inPfp.shw.truth.p.startE;
    ret.shw.truth.p.endE = inPfp.shw.truth.p.endE;
    ret.shw.truth.p.genT = inPfp.shw.truth.p.genT;
    ret.shw.truth.p.startT = inPfp.shw.truth.p.startT;
    ret.shw.truth.p.endT = inPfp.shw.truth.p.endT;
    ret.shw.truth.p.length = inPfp.shw.truth.p.length;
    ret.shw.truth.p.genp.x = inPfp.shw.truth.p.genp.x;
    ret.shw.truth.p.genp.y = inPfp.shw.truth.p.genp.y;
    ret.shw.truth.p.genp.z = inPfp.shw.truth.p.genp.z;
    ret.shw.truth.p.startp.x = inPfp.shw.truth.p.startp.x;
    ret.shw.truth.p.startp.y = inPfp.shw.truth.p.startp.y;
    ret.shw.truth.p.startp.z = inPfp.shw.truth.p.startp.z;
    ret.shw.truth.p.endp.x = inPfp.shw.truth.p.endp.x;
    ret.shw.truth.p.endp.y = inPfp.shw.truth.p.endp.y;
    ret.shw.truth.p.endp.z = inPfp.shw.truth.p.endp.z;
    ret.shw.truth.p.gen.x = inPfp.shw.truth.p.gen.x;
    ret.shw.truth.p.gen.y = inPfp.shw.truth.p.gen.y;
    ret.shw.truth.p.gen.z = inPfp.shw.truth.p.gen.z;
    ret.shw.truth.p.start.x = inPfp.shw.truth.p.start.x;
    ret.shw.truth.p.start.y = inPfp.shw.truth.p.start.y;
    ret.shw.truth.p.start.z = inPfp.shw.truth.p.start.z;
    ret.shw.truth.p.end.x = inPfp.shw.truth.p.end.x;
    ret.shw.truth.p.end.y = inPfp.shw.truth.p.end.y;
    ret.shw.truth.p.end.z = inPfp.shw.truth.p.end.z;
    ret.shw.truth.p.wallin = inPfp.shw.truth.p.wallin;
    ret.shw.truth.p.wallout = inPfp.shw.truth.p.wallout;
    ret.shw.truth.p.cont_tpc = inPfp.shw.truth.p.cont_tpc;
    ret.shw.truth.p.crosses_tpc = inPfp.shw.truth.p.crosses_tpc;
    ret.shw.truth.p.contained = inPfp.shw.truth.p.contained;
    ret.shw.truth.p.pdg = inPfp.shw.truth.p.pdg;
    ret.shw.truth.p.G4ID = inPfp.shw.truth.p.G4ID;
    ret.shw.truth.p.interaction_id = inPfp.shw.truth.p.interaction_id;
    ret.shw.truth.p.cryostat = inPfp.shw.truth.p.cryostat;
    for ( auto const& daughter : inPfp.shw.truth.p.daughters ) ret.shw.truth.p.daughters.push_back( (unsigned)daughter );
    ret.shw.truth.p.parent = inPfp.shw.truth.p.parent;
    ret.shw.truth.p.generator = inPfp.shw.truth.p.generator;
    ret.shw.truth.p.start_process = inPfp.shw.truth.p.start_process;
    ret.shw.truth.p.end_process = inPfp.shw.truth.p.end_process;
    ret.shw.truth.p.gstatus = inPfp.shw.truth.p.gstatus;
    // --> Leave plane, razzle, and selVars info NaN'ed/zero'ed out...

    return;
  }

  void TrackSplitSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    if ( sr->is_clear_cosmic ) return;

    std::vector< caf::SRPFP > retPfps;

    // Vectors for the second tracks...
    std::map< unsigned int, std::vector<std::vector<caf::SRCaloPoint>> > scdryInitIdxToPfpCaloPts;
    std::map< unsigned int, ana::PreservedInitialTrkResults > scdryInitIdxToPreservedInfo;
    unsigned int expectedScdryTrks = 0;

    int pfpIdx = -1;

    // Loop through each pfp and if we should split it, then do the things...
    for ( auto& pfp : sr->reco.pfp ) {
      pfpIdx+=1;
      // Does track cross the cathode or z=0?
      bool crossesCath = false;
      bool crossesMidZ = false;
      if ( std::isnan(pfp.trk.len) || std::isnan(pfp.trk.start.x) ||
           std::isnan(pfp.trk.start.z) || std::isnan(pfp.trk.end.x) || std::isnan(pfp.trk.end.z) ) {
        //caf::SRPFP newPfp;
        //TrackSplitSyst::FillPtrPFP(newPfp, pfp);
        //retPfps.push_back( newPfp );
        continue;
      }
      if ( ((fabs(pfp.trk.start.x)-210.2)*(fabs(pfp.trk.end.x)-210.2)) < 0. )
        crossesCath = true;
      if ( (pfp.trk.start.z * pfp.trk.end.z) < 0. )
        crossesMidZ = true;

      // If both, let's just consider cathode here (TODO: is there a better way?)
      if ( crossesCath && crossesMidZ )
        crossesMidZ = false;
      
      if ( !crossesCath && !crossesMidZ ) {
        //caf::SRPFP newPfp;
        //TrackSplitSyst::FillPtrPFP(newPfp, pfp);
        //retPfps.push_back( newPfp );
        continue;
      }

      // Should we split the track?
      double pToSplit = 0.;
      TVector3 trkAvgDir( pfp.trk.end.x - pfp.trk.start.x, pfp.trk.end.y - pfp.trk.start.y, pfp.trk.end.z - pfp.trk.start.z );
      double nSlopes = 0.;
      if ( crossesCath ) {
        bool isWest = pfp.trk.start.x > 0;
        nSlopes = (isWest ? (210.2-pfp.trk.start.x)/trkAvgDir.X() : (pfp.trk.start.x+210.2)/trkAvgDir.X() );
        double zAtCathode = pfp.trk.start.z + nSlopes*trkAvgDir.Z();
        TVector3 vNormPlane( 1., 0., 0. );
        double ThCathode = trkAvgDir.Angle(vNormPlane);
        ThCathode = std::min(ThCathode, TMath::Pi()-ThCathode);
        double cosThCathode = TMath::Cos( ThCathode );

        const int xBin = splitProb[1][0]->GetXaxis()->FindBin(zAtCathode);
        const int yBin = splitProb[1][0]->GetYaxis()->FindBin(cosThCathode);
        // Had [1][1] here... Believe this is a bug caught 11 March 2024...
        if (xBin == 0 || xBin == splitProb[1][0]->GetXaxis()->GetNbins() + 1 ||
            yBin == 0 || yBin == splitProb[1][0]->GetYaxis()->GetNbins() + 1) {
          //if ( fDebug ) std::cout << "Note that z = " << zAtCathode << " and cosTh = " << cosThCathode << " for cathode." << std::endl;
          //if ( fDebug ) std::cout << "BINS: x=" << xBin << ", y=" << yBin << std::endl;
          //if ( fDebug ) std::cout << "   start: (" << pfp.trk.start.x << ", " << pfp.trk.start.y << ", " << pfp.trk.start.z << ")" << std::endl;
          //if ( fDebug ) std::cout << "   end: (" << pfp.trk.end.x << ", " << pfp.trk.end.y << ", " << pfp.trk.end.z << ")" << std::endl;
          pToSplit = 0;
        }
        else pToSplit = (isWest ? splitProb[1][0]->GetBinContent(xBin,yBin) :
                         splitProb[0][0]->GetBinContent(xBin,yBin) );
      }
      else {
        bool isWest = pfp.trk.start.x > 0;
        bool isNorth = pfp.trk.start.z > 0;
        nSlopes = (!isNorth ? (0.-pfp.trk.start.z)/trkAvgDir.Z() : (pfp.trk.start.z)/trkAvgDir.Z() );
        double xAtMiddleZ = pfp.trk.start.x + nSlopes*trkAvgDir.X();
        TVector3 vNormPlane( 0., 0., 1. );
        double ThMiddleZ = trkAvgDir.Angle(vNormPlane);
        ThMiddleZ = std::min(ThMiddleZ, TMath::Pi()-ThMiddleZ);
        double cosThMiddleZ = TMath::Cos( ThMiddleZ );

        const int xBin = (isWest ? splitProb[1][1]->GetXaxis()->FindBin(xAtMiddleZ) :
                                   splitProb[0][1]->GetXaxis()->FindBin(xAtMiddleZ));
        const int yBin = splitProb[1][1]->GetYaxis()->FindBin(cosThMiddleZ);
        if (xBin == 0 || xBin == splitProb[1][1]->GetXaxis()->GetNbins() + 1 ||
            yBin == 0 || yBin == splitProb[1][1]->GetYaxis()->GetNbins() + 1) {
          //if ( fDebug ) std::cout << "Note that x = " << xAtMiddleZ << " and cosTh = " << cosThMiddleZ << " for z=0 gap." << std::endl;
          //if ( fDebug ) std::cout << "BINS: x=" << xBin << ", y=" << yBin << std::endl;
          //if ( fDebug ) std::cout << "   start: (" << pfp.trk.start.x << ", " << pfp.trk.start.y << ", " << pfp.trk.start.z << ")" << std::endl;
          //if ( fDebug ) std::cout << "   end: (" << pfp.trk.end.x << ", " << pfp.trk.end.y << ", " << pfp.trk.end.z << ")" << std::endl;
          pToSplit = 0;
        }
        else pToSplit = (isWest ? splitProb[1][1]->GetBinContent(xBin,yBin) :
                                  splitProb[0][1]->GetBinContent(xBin,yBin) );
      }
      // DEBUGGING!
      //if ( fDebug ) std::cout << "This track crosses a boundary and its p for splitting is: " << pToSplit << std::endl;
      if ( tRand->Rndm() > pToSplit*sigma ) {
        //if ( fDebug ) std::cout << "... but we did not split it." << std::endl;
        //caf::SRPFP newPfp;
        //TrackSplitSyst::FillPtrPFP(newPfp, pfp);
        //retPfps.push_back( newPfp );
        continue;
      }

      // DEBUGGING!
      //if ( fDebug ) std::cout << "Where? ==> " << (crossesCath ? "Cathode" : "MiddleZ") << std::endl;
      //if ( fDebug ) std::cout << "Trk Start (x,y,z) = (" << pfp.trk.start.x << "," << pfp.trk.start.y << "," << pfp.trk.start.z << ")" << std::endl;
      //if ( fDebug ) std::cout << "Trk End (x,y,z) = (" << pfp.trk.end.x << "," << pfp.trk.end.y << "," << pfp.trk.end.z << ")" << std::endl;

      // DEBUGGING
      //if ( fDebug ) std::cout << "... and we've decided to split the track." << std::endl;

      // Find the calo point at which to break the track on each plane, and save the remaining calo points
      caf::SRCaloPoint splitPoint[3];
      caf::SRCaloPoint firstPoint[3];

      std::vector< std::vector<caf::SRCaloPoint> > savedCaloPoints;

      std::vector< std::vector<caf::SRCaloPoint> > savedScdryCaloPoints;

      ana::PreservedInitialTrkResults PreservedTrkResults;

      for ( unsigned int idxPlane=0; idxPlane < 3; ++idxPlane ){
        savedCaloPoints.push_back({});
        savedScdryCaloPoints.push_back({});
        double minDist = std::numeric_limits<double>::max();
        double largestRR = 0.;
        unsigned int nRRThatAreNaN = 0;
        unsigned int nPoints = 0;

        //unsigned int initialPoints = pfp.trk.calo[idxPlane].points.size();

        for ( auto const& caloPt : pfp.trk.calo[idxPlane].points ) {
          // TODO: REALLY THESE SHOULD GRAB CLOSEST POINT *ON* THE SIDE OF THE SPLIT. THIS IS JUST CLOSEST POINT OVERALL...
          if ( crossesCath ) {
            ////if ( fDebug ) std::cout << caloPt.x << "(" << fabs(caloPt.x)-210.2 << ")";
            if ( fabs(fabs(caloPt.x)-210.2) < minDist ){
              minDist = fabs(fabs(caloPt.x)-210.2);
              splitPoint[idxPlane] = FillCaloPointFrom( caloPt );
              ////if ( fDebug ) std::cout << "[" << splitPoint[idxPlane].x << " " << caloPt.x << "] ";
            }
          }
          else {
            ////if ( fDebug ) std::cout << caloPt.z << "(" << fabs(caloPt.z) << ")";
            if ( fabs(caloPt.z) < minDist ){
              minDist = fabs(caloPt.z);
              splitPoint[idxPlane] = FillCaloPointFrom( caloPt );
              ////if ( fDebug ) std::cout << "[" << splitPoint[idxPlane].z << " " << caloPt.z << "] ";
            }
          }

          if ( !std::isnan(caloPt.rr) && caloPt.rr > largestRR ) {
            largestRR = caloPt.rr;
            firstPoint[idxPlane] = FillCaloPointFrom( caloPt );
          }
          else if ( std::isnan(caloPt.rr) ) {
            nRRThatAreNaN+=1;
          }
        } // Get the split point and first point
        ////if ( fDebug ) std::cout << std::endl;

        for ( auto const& caloPt : pfp.trk.calo[idxPlane].points ) {
          caf::SRCaloPoint savedCaloPoint = FillCaloPointFrom(caloPt);

          if ( caloPt.rr >= splitPoint[idxPlane].rr ) {
            // this is before the split
            savedCaloPoint.rr = caloPt.rr - splitPoint[idxPlane].rr;
            savedCaloPoints[idxPlane].push_back(savedCaloPoint);
            if ( savedCaloPoint.dedx < 1000. ) nPoints+=1;
          }
          else {
            // beyond split point, save to savedScdryCaloPoints
            savedScdryCaloPoints[idxPlane].push_back( savedCaloPoint );
          } 
        } // calo points to save

        // Update the calo info as needed
        pfp.trk.calo[idxPlane].nhit = nPoints;
        pfp.trk.calo[idxPlane].ke = caf::kSignalingNaN;
        pfp.trk.calo[idxPlane].charge = caf::kSignalingNaN;
        pfp.trk.calo[idxPlane].points = savedCaloPoints[idxPlane];

        // The Chi2PID
        TrkChi2Results chi2output = CalculateChi2(pfp.trk.calo[idxPlane]);
        pfp.trk.chi2pid[idxPlane].chi2_proton = chi2output.chi2_proton;
        pfp.trk.chi2pid[idxPlane].chi2_kaon = chi2output.chi2_kaon;
        pfp.trk.chi2pid[idxPlane].chi2_pion = chi2output.chi2_pion;
        pfp.trk.chi2pid[idxPlane].chi2_muon = chi2output.chi2_muon;
        pfp.trk.chi2pid[idxPlane].pida = chi2output.pida;
        pfp.trk.chi2pid[idxPlane].pid_ndof = chi2output.pid_ndof;

        // DEBUGS
        //if ( fDebug ) std::cout << "... PLANE DEBUG INFO: " << idxPlane << std::endl;
        //if ( fDebug ) std::cout << "    initial points = " << initialPoints << std::endl;
        //if ( fDebug ) std::cout << "    --> Post systematic shift" << std::endl;
        //if ( fDebug ) std::cout << "    nRRThatAreNaN = " << nRRThatAreNaN << std::endl;
        //if ( fDebug ) std::cout << "    nhit = " << nPoints << std::endl;
        //if ( fDebug ) std::cout << "    first point (x,y,z) = (" << firstPoint[idxPlane].x << "," << firstPoint[idxPlane].y << "," << firstPoint[idxPlane].z << ")" << std::endl;
        //if ( fDebug ) std::cout << "    first point rr = " << firstPoint[idxPlane].rr << std::endl;
        //if ( fDebug ) std::cout << "    split point (x,y,z) = (" << splitPoint[idxPlane].x << "," << splitPoint[idxPlane].y << "," << splitPoint[idxPlane].z << ")" << std::endl;
        //if ( fDebug ) std::cout << "    split point rr = " << splitPoint[idxPlane].rr << std::endl;
      } // loop planes

      // Save the Preserved Track Info
      PreservedTrkResults.producer = pfp.trk.producer;
      PreservedTrkResults.end_x = pfp.trk.end.x;
      PreservedTrkResults.end_y = pfp.trk.end.y;
      PreservedTrkResults.end_z = pfp.trk.end.z;
      PreservedTrkResults.dir_end_x = pfp.trk.dir_end.x;
      PreservedTrkResults.dir_end_y = pfp.trk.dir_end.y;
      PreservedTrkResults.dir_end_z = pfp.trk.dir_end.z;

      // Update the track info
      pfp.trk.len = ((firstPoint[2].rr-splitPoint[2].rr) + (firstPoint[1].rr-splitPoint[1].rr))/2.;
      pfp.trk.npts = (pfp.trk.calo[2].nhit + pfp.trk.calo[1].nhit + pfp.trk.calo[0].nhit);
      pfp.trk.bestplane = ( pfp.trk.calo[2].nhit >= pfp.trk.calo[1].nhit ? (pfp.trk.calo[2].nhit >= pfp.trk.calo[0].nhit ? caf::Plane_t(2) : caf::Plane_t(0) ) : 
                                                                           (pfp.trk.calo[0].nhit >= pfp.trk.calo[1].nhit ? caf::Plane_t(0) : caf::Plane_t(1) ));

      ////if ( fDebug ) std::cout << "... and its length post-split is now: " << pfp.trk.len << std::endl;

      caf::SRVector3D trkEnd( (splitPoint[2].x + splitPoint[1].x)/2.,
                              (splitPoint[2].y + splitPoint[1].y)/2.,
                              (splitPoint[2].z + splitPoint[1].z)/2. );
      pfp.trk.end = trkEnd;

      caf::SRVector3D trkEndDir( caf::kSignalingNaN, caf::kSignalingNaN, caf::kSignalingNaN );
      pfp.trk.dir_end = trkEndDir;

      TrkMomentumResults p_output = CalculateMomenta((float)pfp.trk.len);
      pfp.trk.rangeP.p_muon = p_output.p_muon;
      pfp.trk.rangeP.p_proton = p_output.p_proton;
      pfp.trk.rangeP.p_pion = p_output.p_pion;

      // NaN some of the remaining variables we don't want the user to use with this shift...
      // TODO: maybe think of rerunning MCS with calo point info?
      caf::SRTrkMCS newMCS;
      pfp.trk.mcsP = newMCS;

      caf::SRTrackScatterClosestApproach newScatter;
      pfp.trk.scatterClosestApproach = newScatter;

      caf::SRTrackStoppingChi2Fit newStopChi2;
      pfp.trk.stoppingChi2Fit = newStopChi2;

      caf::SRTrackDazzle newDazzle;
      pfp.trk.dazzle = newDazzle;

      // Update shower variables... for this one, we ONLY update the shower length to be 90% the track length...
      // -- We may want to update others later...
      pfp.shw.len = 0.9*pfp.trk.len;

      // Save the secondary track info
      unsigned int pfpIdxUI = (unsigned int)pfpIdx;
      scdryInitIdxToPfpCaloPts[ pfpIdxUI ] = savedScdryCaloPoints;
      scdryInitIdxToPreservedInfo[ pfpIdxUI ] = PreservedTrkResults;
      expectedScdryTrks+=1;

      // Fill the new PFP
      caf::SRPFP newPfp;
      TrackSplitSyst::FillPtrPFP(newPfp, pfp);
      retPfps.push_back( newPfp );
    } // loop pfps

    // DEAL WITH SECONDARY TRACKS...
    if ( expectedScdryTrks > 0 )
      //if ( fDebug ) std::cout << "... Dealing with " << expectedScdryTrks << " secondary tracks." << std::endl;
    for ( auto const& [initIdx, caloPtsVec] : scdryInitIdxToPfpCaloPts ) {
      caf::SRPFP newPfp;

      // PFP Elements:
      newPfp.id = -999; // Some default...
      newPfp.parent = sr->reco.pfp[initIdx].id;
      newPfp.parent_is_primary = false;
      newPfp.trackScore = sr->reco.pfp[initIdx].trackScore;
      // --> also the score vars...
      newPfp.pfochar.chgendfrac = sr->reco.pfp[initIdx].pfochar.chgendfrac;
      newPfp.pfochar.chgfracspread = sr->reco.pfp[initIdx].pfochar.chgfracspread;
      newPfp.pfochar.linfitdiff = sr->reco.pfp[initIdx].pfochar.linfitdiff;
      newPfp.pfochar.linfitlen = sr->reco.pfp[initIdx].pfochar.linfitlen;
      newPfp.pfochar.linfitgaplen = sr->reco.pfp[initIdx].pfochar.linfitgaplen;
      newPfp.pfochar.linfitrms = sr->reco.pfp[initIdx].pfochar.linfitrms;
      newPfp.pfochar.openanglediff = sr->reco.pfp[initIdx].pfochar.openanglediff;
      newPfp.pfochar.pca2ratio = sr->reco.pfp[initIdx].pfochar.pca2ratio;
      newPfp.pfochar.pca3ratio = sr->reco.pfp[initIdx].pfochar.pca3ratio;
      newPfp.pfochar.vtxdist = sr->reco.pfp[initIdx].pfochar.vtxdist;
      /////////////////////////////
      newPfp.slcID = sr->reco.pfp[initIdx].slcID;
      newPfp.t0 = (std::isnan(sr->reco.pfp[initIdx].t0) ? caf::kSignalingNaN : (float)sr->reco.pfp[initIdx].t0);

      // Track Elements:
      // Things unchanged from initial track...
      newPfp.trk.producer = scdryInitIdxToPreservedInfo[ initIdx ].producer;

      caf::SRVector3D trkEnd( scdryInitIdxToPreservedInfo[ initIdx ].end_x,
                              scdryInitIdxToPreservedInfo[ initIdx ].end_y,
                              scdryInitIdxToPreservedInfo[ initIdx ].end_z );
      newPfp.trk.end = trkEnd;

      caf::SRVector3D trkEndDir( scdryInitIdxToPreservedInfo[ initIdx ].dir_end_x,
                                 scdryInitIdxToPreservedInfo[ initIdx ].dir_end_y,
                                 scdryInitIdxToPreservedInfo[ initIdx ].dir_end_z );
      newPfp.trk.dir_end = trkEndDir;

      // Things changed from initial track...
      std::vector< std::map<double, caf::SRCaloPoint> > caloPointRRMap;

      for ( unsigned int idxPlane=0; idxPlane < 3; ++idxPlane ){
        caloPointRRMap.push_back( std::map<double, caf::SRCaloPoint>() );
        unsigned int nRRThatAreNaN = 0;
        unsigned int nPoints = 0;

        for ( auto const& caloPt : caloPtsVec[idxPlane] ) {
          if ( !std::isnan(caloPt.rr) ) {
            caloPointRRMap[idxPlane][ caloPt.rr ] = caloPt;
          }
          else if ( std::isnan(caloPt.rr) ) {
            nRRThatAreNaN+=1;
          }
        } // Get the ordering of points to be able to find "first" 5

        for ( auto const& caloPt : caloPtsVec[idxPlane] ) {
          if ( caloPt.dedx < 1000. ) nPoints+=1;
        } // calo points to save

        // Update the calo info as needed
        newPfp.trk.calo[idxPlane].nhit = nPoints;
        newPfp.trk.calo[idxPlane].ke = caf::kSignalingNaN;
        newPfp.trk.calo[idxPlane].charge = caf::kSignalingNaN;
        newPfp.trk.calo[idxPlane].points = caloPtsVec[idxPlane];

        // The Chi2PID
        TrkChi2Results scdryChi2output = CalculateChi2(newPfp.trk.calo[idxPlane]);
        newPfp.trk.chi2pid[idxPlane].chi2_proton = scdryChi2output.chi2_proton;
        newPfp.trk.chi2pid[idxPlane].chi2_kaon = scdryChi2output.chi2_kaon;
        newPfp.trk.chi2pid[idxPlane].chi2_pion = scdryChi2output.chi2_pion;
        newPfp.trk.chi2pid[idxPlane].chi2_muon = scdryChi2output.chi2_muon;
        newPfp.trk.chi2pid[idxPlane].pida = scdryChi2output.pida;
        newPfp.trk.chi2pid[idxPlane].pid_ndof = scdryChi2output.pid_ndof;
      } // loop planes

      //if ( fDebug ) {
      //  std::cout << "Secondary track RRs: ";
      //  for ( auto const& [rrs, pt] : caloPointRRMap[2] ) {
      //    std::cout << rrs << " ";
      //  }
      //  std::cout << std::endl;
      //}

      // Attempt to get the start direction
      std::vector< std::vector<double> > startPoints;
      std::vector< std::vector<double> > dirs;
      std::vector<double> scdryLengths;
      for ( unsigned int idxPlane=1; idxPlane < 3; ++idxPlane ) {
        unsigned int nPoints = 0;
        std::map<double, caf::SRCaloPoint>::iterator it = caloPointRRMap[idxPlane].end();
        if ( caloPointRRMap[idxPlane].size() == 0 ) continue; // if no points then we skip this plane
        --it; // get away from .end()
        scdryLengths.push_back( it->second.rr );
        std::vector<double> startPoint = { it->second.x, it->second.y, it->second.z };
        startPoints.push_back(startPoint);
        while ( it!=caloPointRRMap[idxPlane].begin() ) {
          double xyz_start[3] = { it->second.x, it->second.y, it->second.z };
          //std::cout << "xyz: " << xyz_start[0] << " " << xyz_start[1] << " " << xyz_start[2] << std::endl;
          --it;
          double xyz_end[3] = { it->second.x, it->second.y, it->second.z };
          nPoints+=1;
          std::vector<double> dir = { xyz_end[0]-xyz_start[0], xyz_end[1]-xyz_start[1], xyz_end[2]-xyz_start[2] };
          dirs.push_back(dir);
          if ( nPoints==5 ) break;
        }
      }

      std::vector< double > ave_dir = {0., 0., 0.};
      for ( auto const& dir: dirs ) {
        TVector3 dir_tv3( dir[0], dir[1], dir[2] );
        dir_tv3 = dir_tv3.Unit();
        ave_dir[0] += dir_tv3.X();
        ave_dir[1] += dir_tv3.Y();
        ave_dir[2] += dir_tv3.Z();
      }
      if ( dirs.size() != 0 ) {
        ave_dir[0] /= double(dirs.size());
        ave_dir[1] /= double(dirs.size());
        ave_dir[2] /= double(dirs.size());
      }

      TVector3 ave_dir_tv3( ave_dir[0], ave_dir[1], ave_dir[2] );
      if ( dirs.size() != 0 ) ave_dir_tv3 = ave_dir_tv3.Unit();

      caf::SRVector3D trkStartDir( ave_dir_tv3.X(),
                                   ave_dir_tv3.Y(),
                                   ave_dir_tv3.Z() );
      newPfp.trk.dir = trkStartDir;

      // Length
      double scdryLength = 0.;
      for ( auto const& scdryLenVal : scdryLengths ) scdryLength+=scdryLenVal;
      if( scdryLengths.size() != 0 ) scdryLength /= double(scdryLengths.size());
      newPfp.trk.len = scdryLength;

      // Start point
      double startPtX = 0.;
      double startPtY = 0.;
      double startPtZ = 0.;
      for ( auto const& startPoint : startPoints ) {
        startPtX+=startPoint[0];
        startPtY+=startPoint[1];
        startPtZ+=startPoint[2];
      }
      if ( startPoints.size() != 0 ) {
        startPtX /= double(startPoints.size());
        startPtY /= double(startPoints.size());
        startPtZ /= double(startPoints.size());
        caf::SRVector3D trkStartLoc( startPtX,
                                     startPtY,
                                     startPtZ );
        newPfp.trk.start = trkStartLoc;
      }
      else {
        caf::SRVector3D trkStartLoc( caf::kSignalingNaN,
                                     caf::kSignalingNaN,
                                     caf::kSignalingNaN );
        newPfp.trk.start = trkStartLoc;
      }

      // N points & best plane
      newPfp.trk.npts = (newPfp.trk.calo[2].nhit + newPfp.trk.calo[1].nhit + newPfp.trk.calo[0].nhit);
      newPfp.trk.bestplane = ( newPfp.trk.calo[2].nhit >= newPfp.trk.calo[1].nhit ? (newPfp.trk.calo[2].nhit >= newPfp.trk.calo[0].nhit ? caf::Plane_t(2) : caf::Plane_t(0) ) : 
                                                                                    (newPfp.trk.calo[0].nhit >= newPfp.trk.calo[1].nhit ? caf::Plane_t(0) : caf::Plane_t(1) ));

      TrkMomentumResults p_output = CalculateMomenta(newPfp.trk.len);
      newPfp.trk.rangeP.p_muon = p_output.p_muon;
      newPfp.trk.rangeP.p_proton = p_output.p_proton;
      newPfp.trk.rangeP.p_pion = p_output.p_pion;

      // Shower Elements:
      // for our shower cut we check shw.start.x, shw.len, shw.conversion_gap
      // Update length as above
      newPfp.shw.len = 0.9*newPfp.trk.len;
      // Make start X, Y, Z be the same as the track start
      newPfp.shw.start.x = !std::isnan(newPfp.trk.start.x) ? newPfp.trk.start.x : caf::kSignalingNaN;
      newPfp.shw.start.y = !std::isnan(newPfp.trk.start.y) ? newPfp.trk.start.y : caf::kSignalingNaN;
      newPfp.shw.start.z = !std::isnan(newPfp.trk.start.z) ? newPfp.trk.start.z : caf::kSignalingNaN;
      // Make conversion gap be the distance between the slice vertex and this point...
      if ( std::isnan(newPfp.shw.start.x) || std::isnan(newPfp.shw.start.y) || std::isnan(newPfp.shw.start.z) ) {
        newPfp.shw.conversion_gap = caf::kSignalingNaN;
      }
      else {
        newPfp.shw.conversion_gap = std::hypot( newPfp.shw.start.x - sr->vertex.x,
                                                newPfp.shw.start.y - sr->vertex.y,
                                                newPfp.shw.start.z - sr->vertex.z );
      }

      // ... and for now leave EVERYTHING else alone.
      // Buyer beware of this. If other groups/analyses want to use this, we may need to rethink.
      retPfps.push_back( newPfp ); // ORIGINAL
      //sr->reco.pfp.push_back( newPfp ); // TESTING!!
    } // secondary tracks

    //     std::map< unsigned int, ana::PreservedInitialTrkResults > scdryInitIdxToPreservedInfo;
    if ( scdryInitIdxToPreservedInfo.size()==0 ) {
      // No split tracks
      return;
    } // return NO SPLIT
    else {
      // We have the tracks that were split put into new retPfps. We need to make the OTHERS into caf::SRPFP objects.
      for ( unsigned int idxPfp = 0; idxPfp < sr->reco.pfp.size(); ++idxPfp ) {
        bool isThisTheSplitTrack = false;
        for ( auto const& [idxSplit, preservedInfo] : scdryInitIdxToPreservedInfo ) {
          if ( idxPfp == idxSplit ) {
            isThisTheSplitTrack = true;
            break;
          }
        } // loop split tracks
        if ( isThisTheSplitTrack ) continue;

        // Fill a new PFP for this PFP
        caf::SRPFP newPfp;
        TrackSplitSyst::FillPtrPFP(newPfp, sr->reco.pfp[idxPfp]);
        retPfps.push_back( newPfp );
      } // loop pfps

      // Update the slice reco info with these SRPFPs
      sr->reco.npfp = retPfps.size();
      sr->reco.pfp = retPfps;

      return;
    } // return WITH SPLIT

    // Update the slice reco info with these SRPFPs
    //sr->reco.npfp = retPfps.size();
    //sr->reco.pfp = retPfps;


  } // Shift

  void TrackSplitSyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

  }

  const TrackSplitSyst kTrackSplittingSyst("TrackSplittingSyst", "Split tracks systematic");
  const TrackSplitSyst kTrackSplittingSystDebug("TrackSplittingSyst", "Split tracks systematic", true);



  // Class for checking the approximations/calculations used in the systematic
  void TrackSplitSystCheck::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    if ( sr->is_clear_cosmic ) return;
    if ( sr->truth.index < 0 ) return;

    // Let's just do this for PFPs with tracks that are longer than 50cm and match to true muons?
    for ( unsigned int idxPfp = 0; idxPfp < sr->reco.npfp; ++idxPfp ) {
      if ( std::isnan(sr->reco.pfp[idxPfp].trk.truth.p.pdg) || std::isinf(sr->reco.pfp[idxPfp].trk.truth.p.pdg) ||
           std::isnan(sr->reco.pfp[idxPfp].trk.start.x) || std::isinf(sr->reco.pfp[idxPfp].trk.start.x) ||
           std::isnan(sr->reco.pfp[idxPfp].trk.len) || std::isinf(sr->reco.pfp[idxPfp].trk.len) ||
           abs(sr->reco.pfp[idxPfp].trk.truth.p.pdg) != 13 || sr->reco.pfp[idxPfp].trk.len < 50. || sr->reco.pfp[idxPfp].trk.calo[2].nhit < 10 || 
           std::isnan(sr->reco.pfp[idxPfp].shw.conversion_gap) || std::isinf(sr->reco.pfp[idxPfp].shw.conversion_gap) ) {
        continue;
      }

      caf::SRCaloPoint lastPoint[3];
      caf::SRCaloPoint firstPoint[3];
      std::vector< std::map<double, caf::SRCaloPoint> > caloPointRRMap;

      for ( unsigned int idxPlane=0; idxPlane < 3; ++idxPlane ){
        double minDist = std::numeric_limits<double>::max();
        double largestRR = 0.;
        unsigned int nRRThatAreNaN = 0;
        caloPointRRMap.push_back( std::map<double, caf::SRCaloPoint>() );

        for ( auto const& caloPt : sr->reco.pfp[idxPfp].trk.calo[idxPlane].points ) {

          // First and last points for start/end position and length and such
          if ( !std::isnan(caloPt.rr) ) {
            caloPointRRMap[idxPlane][ caloPt.rr ] = FillCaloPointFrom(caloPt);

            if ( caloPt.rr < minDist ) {
              minDist = caloPt.rr;
              lastPoint[idxPlane] = FillCaloPointFrom( caloPt );
            }
            if ( caloPt.rr > largestRR ) {
              largestRR = caloPt.rr;
              firstPoint[idxPlane] = FillCaloPointFrom( caloPt );
            }
          }
          else {
            nRRThatAreNaN+=1;
          }
        } // loop calo points

        // Do the chi2PID calculation --> didn't change anything here so hopefully see no change...
        TrkChi2Results chi2output = CalculateChi2(sr->reco.pfp[idxPfp].trk.calo[idxPlane]);
        sr->reco.pfp[idxPfp].trk.chi2pid[idxPlane].chi2_proton = chi2output.chi2_proton;
        sr->reco.pfp[idxPfp].trk.chi2pid[idxPlane].chi2_kaon = chi2output.chi2_kaon;
        sr->reco.pfp[idxPfp].trk.chi2pid[idxPlane].chi2_pion = chi2output.chi2_pion;
        sr->reco.pfp[idxPfp].trk.chi2pid[idxPlane].chi2_muon = chi2output.chi2_muon;
        sr->reco.pfp[idxPfp].trk.chi2pid[idxPlane].pida = chi2output.pida;
        sr->reco.pfp[idxPfp].trk.chi2pid[idxPlane].pid_ndof = chi2output.pid_ndof;
      } // loop plane

      // Things to get "edited" versions of
      // NOTE: temporarily commenting out the length updating so that I can get the same tracks before and after editing by picking on length...
      //// sr->reco.pfp[idxPfp].trk.len = ((firstPoint[2].rr-lastPoint[2].rr) + (firstPoint[1].rr-lastPoint[1].rr))/2.; // TEMPORARILY COMMENTED OUT
      // NOTE: we could also think of trying the version used for the 2nd track...
      sr->reco.pfp[idxPfp].trk.npts = (sr->reco.pfp[idxPfp].trk.calo[2].nhit + sr->reco.pfp[idxPfp].trk.calo[1].nhit + sr->reco.pfp[idxPfp].trk.calo[0].nhit);
      sr->reco.pfp[idxPfp].trk.bestplane = ( sr->reco.pfp[idxPfp].trk.calo[2].nhit >= sr->reco.pfp[idxPfp].trk.calo[1].nhit ? (sr->reco.pfp[idxPfp].trk.calo[2].nhit >= sr->reco.pfp[idxPfp].trk.calo[0].nhit ? caf::Plane_t(2) : caf::Plane_t(0) ) :
                                                                                                                              (sr->reco.pfp[idxPfp].trk.calo[0].nhit >= sr->reco.pfp[idxPfp].trk.calo[1].nhit ? caf::Plane_t(0) : caf::Plane_t(1) ));

      // -- Beginning of track
      std::vector< std::vector<double> > startPoints;
      std::vector< std::vector<double> > dirs;
      std::vector<double> scdryLengths;
      for ( unsigned int idxPlane=0; idxPlane < 3; ++idxPlane ) {
        unsigned int nPoints = 0;
        std::map<double, caf::SRCaloPoint>::iterator it = caloPointRRMap[idxPlane].end();
        if ( caloPointRRMap[idxPlane].size() == 0 ) continue; // if no points then we skip this plane
        --it; // get away from .end()
        scdryLengths.push_back( it->second.rr );
        std::vector<double> startPoint = { it->second.x, it->second.y, it->second.z };
        startPoints.push_back(startPoint);
	if ( idxPlane==0 ) continue;
        while ( it!=caloPointRRMap[idxPlane].begin() ) {
          double xyz_start[3] = { it->second.x, it->second.y, it->second.z };
          //std::cout << "xyz: " << xyz_start[0] << " " << xyz_start[1] << " " << xyz_start[2] << std::endl;
          --it;
          double xyz_end[3] = { it->second.x, it->second.y, it->second.z };
          nPoints+=1;
          std::vector<double> dir = { xyz_end[0]-xyz_start[0], xyz_end[1]-xyz_start[1], xyz_end[2]-xyz_start[2] };
          dirs.push_back(dir);
          if ( nPoints==5 ) break;
        }
      }

      std::vector< double > ave_dir = {0., 0., 0.};
      for ( auto const& dir: dirs ) {
        TVector3 dir_tv3( dir[0], dir[1], dir[2] );
        dir_tv3 = dir_tv3.Unit();
        ave_dir[0] += dir_tv3.X();
        ave_dir[1] += dir_tv3.Y();
        ave_dir[2] += dir_tv3.Z();
      }
      if ( dirs.size() != 0 ) {
        ave_dir[0] /= double(dirs.size());
        ave_dir[1] /= double(dirs.size());
        ave_dir[2] /= double(dirs.size());
      }

      TVector3 ave_dir_tv3( ave_dir[0], ave_dir[1], ave_dir[2] );
      if ( dirs.size() != 0 ) ave_dir_tv3 = ave_dir_tv3.Unit();

      caf::SRVector3D trkStartDir( ave_dir_tv3.X(),
                                   ave_dir_tv3.Y(),
                                   ave_dir_tv3.Z() );
      sr->reco.pfp[idxPfp].trk.dir = trkStartDir;

      // Start point
      /*
      double startPtX = 0.;
      double startPtY = 0.;
      double startPtZ = 0.;
      for ( auto const& startPoint : startPoints ) {
        startPtX+=startPoint[0];
        startPtY+=startPoint[1];
        startPtZ+=startPoint[2];
      }
      if ( startPoints.size() != 0 ) {
        startPtX /= double(startPoints.size());
        startPtY /= double(startPoints.size());
        startPtZ /= double(startPoints.size());
        caf::SRVector3D trkStartLoc( startPtX,
                                     startPtY,
                                     startPtZ );
        sr->reco.pfp[idxPfp].trk.start = trkStartLoc;
      }
      else {
        // caf::SRVector3D trkStartLoc( caf::kSignalingNaN,
        //                             caf::kSignalingNaN,
        //                             caf::kSignalingNaN );
        // FOR PURPOSES OF THIS MODULE, LET'S MAKE THIS -9999,-9999,-9999 so the distance is very large but not nans...
        caf::SRVector3D trkStartLoc( -9999,-9999,-9999 );
        sr->reco.pfp[idxPfp].trk.start = trkStartLoc;
      }
      */
      // Instead of averaging the first point on each of ind2 and collection, pick the one with the largest RR and use that directly...
      double startRR = 0.;
      double startPtX = -9999.;
      double startPtY = -9999.;
      double startPtZ = -9999.;
      for ( unsigned int idxPt = 0 ; idxPt < startPoints.size(); ++idxPt ) {
        if ( scdryLengths[idxPt] > startRR ) {
          startRR = scdryLengths[idxPt];
          startPtX = startPoints[idxPt][0];
          startPtY = startPoints[idxPt][1];
          startPtZ = startPoints[idxPt][2];
        }
      }
      caf::SRVector3D trkStartLoc( startPtX,
				   startPtY,
				   startPtZ );
      sr->reco.pfp[idxPfp].trk.start = trkStartLoc;

      // -- End of track
      //caf::SRVector3D trkEnd( (lastPoint[2].x + lastPoint[1].x)/2.,
      //                        (lastPoint[2].y + lastPoint[1].y)/2.,
      //                        (lastPoint[2].z + lastPoint[1].z)/2. );
      //sr->reco.pfp[idxPfp].trk.end = trkEnd;
      // and TRY SIMILAR WITH THE LAST POINT...
      double endRR = std::numeric_limits<double>::max();
      double endPtX = -9999.;
      double endPtY = -9999.;
      double endPtZ = -9999.;
      for ( unsigned int idxPt = 0; idxPt < 3; ++idxPt ) {
      	if ( std::isnan(lastPoint[idxPt].rr) || std::isinf(lastPoint[idxPt].rr) ) continue;
        if ( lastPoint[idxPt].rr < endRR ) {
          endRR = lastPoint[idxPt].rr;
          endPtX = lastPoint[idxPt].x;
          endPtY = lastPoint[idxPt].y;
          endPtZ = lastPoint[idxPt].z;
        }
      }
      caf::SRVector3D trkEnd( endPtX, endPtY, endPtZ );
      sr->reco.pfp[idxPfp].trk.end = trkEnd;
      // -- But actually average over:

      // -- Other track and shower quantities
      TrkMomentumResults p_output = CalculateMomenta(sr->reco.pfp[idxPfp].trk.len);
      sr->reco.pfp[idxPfp].trk.rangeP.p_muon = p_output.p_muon;
      sr->reco.pfp[idxPfp].trk.rangeP.p_proton = p_output.p_proton;
      sr->reco.pfp[idxPfp].trk.rangeP.p_pion = p_output.p_pion;

      // Shower Elements:
      // for our shower cut we check shw.start.x, shw.len, shw.conversion_gap
      // Update length as above
      // *NOTE:* these I do expect to see some differences... but, not sure what else to do at moment without being super complicated...
      ///// sr->reco.pfp[idxPfp].shw.len = 0.9*sr->reco.pfp[idxPfp].trk.len; // TEMPORARILY COMMENTED OUT
      // Make start X, Y, Z be the same as the track start
      sr->reco.pfp[idxPfp].shw.start.x = !std::isnan(sr->reco.pfp[idxPfp].trk.start.x) ? (float)sr->reco.pfp[idxPfp].trk.start.x : caf::kSignalingNaN;
      sr->reco.pfp[idxPfp].shw.start.y = !std::isnan(sr->reco.pfp[idxPfp].trk.start.y) ? (float)sr->reco.pfp[idxPfp].trk.start.y : caf::kSignalingNaN;
      sr->reco.pfp[idxPfp].shw.start.z = !std::isnan(sr->reco.pfp[idxPfp].trk.start.z) ? (float)sr->reco.pfp[idxPfp].trk.start.z : caf::kSignalingNaN;
      // Make conversion gap be the distance between the slice vertex and this point...
      if ( std::isnan(sr->reco.pfp[idxPfp].shw.start.x) || std::isnan(sr->reco.pfp[idxPfp].shw.start.y) || std::isnan(sr->reco.pfp[idxPfp].shw.start.z) ) {
        // FOR PURPOSES OF THIS MODULE, LET'S MAKE THIS -9999 so the distance is very large but not nans...
        //sr->reco.pfp[idxPfp].shw.conversion_gap = caf::kSignalingNaN;
        sr->reco.pfp[idxPfp].shw.conversion_gap = -9999;
      }
      else {
        sr->reco.pfp[idxPfp].shw.conversion_gap = std::hypot( sr->reco.pfp[idxPfp].shw.start.x - sr->vertex.x,
                                                              sr->reco.pfp[idxPfp].shw.start.y - sr->vertex.y,
                                                              sr->reco.pfp[idxPfp].shw.start.z - sr->vertex.z );
      }

    } // loop pfps
  } // Shift

  void TrackSplitSystCheck::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

  }

  const TrackSplitSystCheck kTrackSplittingSystCheck("TrackSplittingSystCheck", "Split tracks systematic testing");


} // end namespace ana
