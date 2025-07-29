#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ana {
  NuMIFluxSyst::NuMIFluxSyst(const std::string& dir,
                             const std::string& prefix,
                             const std::string& name,
                             const std::string& fluxFile /* = "" */)
    : ISyst("numi_" + name, "NuMI flux: " + name)
    , fHistName(dir + "/" + prefix + name)
    , fName(name)
    , fScale()
  {
    if (!fluxFile.empty()) { fFluxFilePath = fluxFile; }
    else {
      const char* sbndata = std::getenv("SBNDATA_DIR");
      if (!sbndata) {
        std::cout << "NuMIPpfxFluxWeight: $SBNDATA_DIR environment variable not set. Please setup "
                     "the sbndata product."
                  << std::endl;
        std::abort();
      }

      fFluxFilePath = std::string(sbndata) +
                     "beamData/NuMIdata/2025-04-08_out_450.37_7991.98_79512.66.root";
    }
  }

  //----------------------------------------------------------------------
  NuMIFluxSyst::~NuMIFluxSyst()
  {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          delete fScale[i][j][k];
  }

  //----------------------------------------------------------------------
  void NuMIFluxSyst::Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const
  {
    return this->Shift(sigma, &slc->truth, weight);
  }

  //----------------------------------------------------------------------
  void NuMIFluxSyst::Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const
  {
    if(nu->index < 0) return; // No neutrino

    // TODO switch to regular detector flag as soon as it's filled reliably
    if(nu->baseline < 500){
      std::cout << "Using NuMIFluxSyst on what appears not to be an icarus file (baseline = " << nu->baseline << ")!" << std::endl;
      abort();
    }

    if (!fScale[0][0][0]) {
      TFile f(fFluxFilePath.c_str());
      if (f.IsZombie()) {
        std::cout << "NuMIFluxSysts: Failed to open " << fFluxFilePath << std::endl;
        abort();
      }

      for (int hcIdx : {0, 1}) {
        for (int flavIdx : {0, 1}) {
          for (int signIdx : {0, 1}) {
            std::string hName = fHistName;
            if (hcIdx == 0)
              hName += "_fhc";
            else
              hName += "_rhc";
            if (flavIdx == 0)
              hName += "_nue";
            else
              hName += "_numu";
            if (signIdx == 1) hName += "bar";

            TH1* h = (TH1*)f.Get(hName.c_str());
            if (!h) {
              std::cout << "NuMIFluxSysts: failed to find " << hName << " in " << f.GetName()
                        << std::endl;
              abort();
            }
            h = (TH1*)h->Clone(UniqueName().c_str());
            h->SetDirectory(0);

            fScale[hcIdx][flavIdx][signIdx] = h;
          }
        }
      }
    } // end if

    const int flav = abs(nu->initpdg);
    if(flav != 12 && flav != 14) return;

    const int hcIdx = 0; // TODO TODO always FHC
    const int flavIdx = (flav == 12) ? 0 : 1;
    const int signIdx = (nu->initpdg > 0) ? 0 : 1;

    TH1* h = fScale[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(nu->E);

    if (bin == 0 || bin == h->GetNbinsX() + 1) return;

    // BH -- attempt to not change any weight if the value is inf or nan, just continue with weight as-is...
    if (std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin)))
      weight *= 1;
    else
      weight *= 1 + sigma * h->GetBinContent(bin);
  }

  //----------------------------------------------------------------------
  const NuMIFluxSyst* GetNuMIFluxSyst(const std::string& dir,
                                      const std::string& prefix,
                                      const std::string& name,
                                      const std::string& fluxFile /* = "" */)
  {
    // Make sure we always give the same one back
    static std::map<std::string, const NuMIFluxSyst*> cache;

    const std::string key = dir + "/" + prefix + name;
    if (cache.count(key) == 0) cache[key] = new NuMIFluxSyst(dir, prefix, name, fluxFile);

    return cache[key];
  }

  const ISyst* GetNuMIFluxStatisticalUncertainty()
  {
    return GetNuMIFluxSyst("statistical_uncertainties", "h", "stat");
  }

  //----------------------------------------------------------------------
  std::vector<const ISyst*> GetNuMIHadronProductionFluxSysts()
  {
    /* Tony W: These should not be used directly for syst shifts; prefer to use the PCA form.
     * The HP fractional uncertainties are only included for completeness.
     */
    const std::vector<std::string> syst_names = {
      "att", "mesinc", "nCpi", "nuAlFe", "nua", "others", "pCnu", "pCpi"};

    std::vector<const ISyst*> ret;
    for (const auto& name : syst_names)
      ret.push_back(
        GetNuMIFluxSyst("fractional_uncertainties/hadron/" + name, "hfrac_hadron_", name));
    return ret;
  }

  //----------------------------------------------------------------------
  // This is valid when using 2023-07-31_out_450.37_7991.98_79512.66_QEL11.root
  // If you are using 2025-04-08_out_450.37_7991.98_79512.66.root, use GetNuMIBeamShiftSysts()
  std::vector<const ISyst*> GetNuMIBeamlineFluxSysts()
  {

    std::cout << "[GetNuMIBeamlineFluxSysts] This is valid for when using 2023-07-31_out_450.37_7991.98_79512.66_QEL11.root" << std::endl;
    std::cout << "[GetNuMIBeamlineFluxSysts] If you are using 2025-04-08_out_450.37_7991.98_79512.66.root, use GetNuMIBeamShiftSysts()" << std::endl;

    const std::vector<std::string> syst_names = {"beam_div",
                                                 "beam_power",
                                                 "beam_shift_y_minus",
                                                 "beam_shift_y_plus",
                                                 "beam_shift_x",
                                                 "beam_spot",
                                                 "horn1_x",
                                                 "horn1_y",
                                                 "horn_current_plus",
                                                 "water_layer"};

    // TODO) USE OLD FILE FOR BEAMLINE UNTIL NEW ONE ARRIVES
    const char* sbndata = std::getenv("SBNDATA_DIR");
    std::string fFluxFilePath = std::string(sbndata) +
                     "beamData/NuMIdata/2023-07-31_out_450.37_7991.98_79512.66_QEL11.root";

    std::vector<const ISyst*> ret;
    for (const auto& name : syst_names)
      ret.push_back(GetNuMIFluxSyst("beam_systematic_shifts", "hsyst_beam_", name, fFluxFilePath));
    return ret;
  }

  //----------------------------------------------------------------------
  std::vector<const ISyst*> GetNuMIPCAFluxSysts(unsigned int Npcs)
  {
    std::vector<const ISyst*> ret;
    for (unsigned int i = 0; i < Npcs; ++i) {
      ret.push_back(GetNuMIFluxSyst("pca/principal_components", "h", "pc_" + std::to_string(i)));
    }
    return ret;
  }

  //----------------------------------------------------------------------
  std::vector<const ISyst*> GetAllNuMIFluxSysts(unsigned int Npcs)
  {
    std::vector<const ISyst*> ret = GetNuMIBeamShiftSysts();
    std::vector<const ISyst*> add = GetNuMIPCAFluxSysts(Npcs);
    const ISyst* stat = GetNuMIFluxStatisticalUncertainty();
    ret.insert(ret.end(), add.begin(), add.end());
    ret.push_back(stat);
    return ret;
  }

  //----------------------------------------------------------------------
  // Flat X percent uncertainty on the G3Chase correction (i.e., the missing steel block)

  NuMIBeamG3ChaseSyst::NuMIBeamG3ChaseSyst(const std::string& name, const std::string& latexName):
    ISyst(name, latexName)
  {

  }

  void NuMIBeamG3ChaseSyst::Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const
  {
    this->Shift(sigma, &sr->truth, weight);
  }

  void NuMIBeamG3ChaseSyst::Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const {

    if(sr->index < 0) return; // No neutrino

    static double G3ChaseOneSigError = 0.086; // 8.6percent
    double this_rw = 1. + sigma * G3ChaseOneSigError;
    weight *= this_rw;

  }

  //----------------------------------------------------------------------
  // NuMI Flux Systematic for BeamShift specifically, combining two hists from the file

  NuMIBeamShiftSyst::NuMIBeamShiftSyst(const std::string& name,
                                       const std::string& fluxFile)
    : ISyst("numi_" + name, "NuMI flux, beam focusing: " + name)
    , fName(name)
    , fScale()
  {
    if (!fluxFile.empty()) { fFluxFilePath = fluxFile; }
    else {
      const char* sbndata = std::getenv("SBNDATA_DIR");
      if (!sbndata) {
        std::cout << "NuMIBeamShiftSyst: $SBNDATA_DIR environment variable not set. Please setup "
                     "the sbndata product."
                  << std::endl;
        std::abort();
      }

      // TODO) UPDATE TO NEWER FILE ONCE BEAM SHIFT SYSTS ARRIVES

      fFluxFilePath = std::string(sbndata) +
                      "beamData/NuMIdata/2025-04-08_out_450.37_7991.98_79512.66.root";
    }
  }

  //----------------------------------------------------------------------
  NuMIBeamShiftSyst::~NuMIBeamShiftSyst()
  {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k)
          for (int l = 0; l < 2; ++l)
            delete fScale[i][j][k][l];
  }

  //----------------------------------------------------------------------
  void NuMIBeamShiftSyst::Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const
  {
    return this->Shift(sigma, &slc->truth, weight);
  }

  //----------------------------------------------------------------------
  void NuMIBeamShiftSyst::Shift(double sigma, caf::SRTrueInteractionProxy* nu, double& weight) const
  {
    if(nu->index < 0) return; // No neutrino

    // TODO switch to regular detector flag as soon as it's filled reliably
    if(nu->baseline < 500){
      std::cout << "Using NuMIFluxSyst on what appears not to be an icarus file (baseline = " << nu->baseline << ")!" << std::endl;
      abort();
    }

    if (!fScale[0][0][0][0]) {
      TFile f(fFluxFilePath.c_str());
      if (f.IsZombie()) {
        std::cout << "NuMIBeamShiftSyst: Failed to open " << fFluxFilePath << std::endl;
        abort();
      }

      for ( int sigmaIdx : {0, 1} ) {

        std::string ShiftHistName = "hsyst_beam";

        if(fName=="HornCurr"){
          if(sigmaIdx==0) ShiftHistName += "_Horn_p2kA";
          else            ShiftHistName += "_Horn_m2kA";
        }
        else if(fName=="Horn1_x"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"p3mm";
          else            ShiftHistName += "_"+fName+"_"+"m3mm";
        }
        else if(fName=="Horn1_y"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"p3mm";
          else            ShiftHistName += "_"+fName+"_"+"m3mm";
        }
        else if(fName=="Beam_spot"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"1_7mm";
          else            ShiftHistName += "_"+fName+"_"+"1_3mm";
        }
        else if(fName=="Horn2_x"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"p3mm";
          else            ShiftHistName += "_"+fName+"_"+"m3mm";
        }
        else if(fName=="Horn2_y"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"p3mm";
          else            ShiftHistName += "_"+fName+"_"+"m3mm";
        }
        else if(fName=="Horns_water"){
          if(sigmaIdx==0) ShiftHistName += "_Horns_2mm_water";
          else            ShiftHistName += "_Horns_0mm_water";
        }
        else if(fName=="Beam_shift_x"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"p1mm";
          else            ShiftHistName += "_"+fName+"_"+"m1mm";
        }
        else if(fName=="Beam_shift_y"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"p1mm";
          else            ShiftHistName += "_"+fName+"_"+"m1mm";
        }
        else if(fName=="Target_z"){
          if(sigmaIdx==0) ShiftHistName += "_"+fName+"_"+"p7mm";
          else            ShiftHistName += "_"+fName+"_"+"m7mm";
        }
        else{
            std::cout << "[NuMIBeamShiftSyst::Shift] " << fName << " is not supported " << std::endl;
            abort();
        }

        for (int hcIdx : {0, 1}) {
          for (int flavIdx : {0, 1}) {
            for (int signIdx : {0, 1}) {

              std::string hName = "";
              if(hcIdx==0) hName = "beam_focusing_uncertainties/fhc/"+ShiftHistName+"_fhc";
              else         hName = "beam_focusing_uncertainties/rhc/"+ShiftHistName+"_rhc";

              if(flavIdx==0) hName += "_nue";
              else           hName += "_numu";

              if (signIdx == 1) hName += "bar";

              // TODO
              // 05/28/2025 JK) Type in current file
              if(fName=="Horn1_x" && sigmaIdx==1){
                hName = "";
                std::string tmp_ShiftHistName = "hsyst_beam_Horm1_x_m3mm";
                if(hcIdx==0) hName = "beam_focusing_uncertainties/fhc/"+tmp_ShiftHistName+"_fhc";
                else         hName = "beam_focusing_uncertainties/rhc/"+tmp_ShiftHistName+"_rhc";

                if(flavIdx==0) hName += "_nue";
                else           hName += "_numu";

                if (signIdx == 1) hName += "bar";                
              }

              TH1* h = (TH1*)f.Get(hName.c_str());
              if (!h) {
                std::cout << "NuMIBeamShiftSyst: failed to find " << hName << " in " << f.GetName()
                          << std::endl;
                abort();
              }
              h = (TH1*)h->Clone(UniqueName().c_str());
              h->SetDirectory(0);

              fScale[sigmaIdx][hcIdx][flavIdx][signIdx] = h;

            } // sign idx
          } // flavor idx
        } // hc idx
      } // sigma idx
    } // end if

    const int flav = abs(nu->initpdg);
    if(flav != 12 && flav != 14) return;

    const int hcIdx = 0; // TODO TODO always FHC
    const int flavIdx = (flav == 12) ? 0 : 1;
    const int signIdx = (nu->initpdg > 0) ? 0 : 1;

    const int sigmaIdx = (sigma >= 0) ? 0 : 1;
    sigma = fabs(sigma);

    TH1* h = fScale[sigmaIdx][hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(nu->E);

    if (bin == 0 || bin == h->GetNbinsX() + 1) return;

    // BH -- attempt to not change any weight if the value is inf or nan, just continue with weight as-is...
    if (std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin)))
      weight *= 1;
    else
      weight *= 1 + sigma * h->GetBinContent(bin);
  }

  const NuMIBeamShiftSyst* GetNuMIBeamShiftSyst(const std::string& name,
                                                const std::string& fluxFile)
  {
    // Make sure we always give the same one back
    static std::map<std::string, const NuMIBeamShiftSyst*> cache;

    const std::string key = name;
    if (cache.count(key) == 0) cache[key] = new NuMIBeamShiftSyst(name, fluxFile);

    return cache[key];
  }

  std::vector<const ISyst*> GetNuMIBeamShiftSysts(){

    const std::vector<std::string> syst_names = {
      "HornCurr",
      "Horn1_x",
      "Horn1_y",
      "Beam_spot",
      "Horn2_x",
      "Horn2_y",
      "Horns_water",
      "Beam_shift_x",
      "Beam_shift_y",
      "Target_z",
    };

    const char* sbndata = std::getenv("SBNDATA_DIR");
    std::string fFluxFilePath = std::string(sbndata) +
                     "beamData/NuMIdata/2025-04-08_out_450.37_7991.98_79512.66.root";

    std::vector<const ISyst*> ret;
    for (const auto& name : syst_names){
      ret.push_back(GetNuMIBeamShiftSyst(name, fFluxFilePath));
    }

    return ret;

  }

} // namespace ana
