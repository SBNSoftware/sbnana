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
                     "beamData/NuMIdata/2023-07-31_out_450.37_7991.98_79512.66_QEL11.root";
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
    if (slc->truth.index < 0) return; // No neutrino

    // TODO switch to regular detector flag as soon as it's filled reliably
    if (slc->truth.baseline < 500) {
      std::cout << "Using NuMIFluxSyst on what appears not to be an icarus file (baseline = "
                << slc->truth.baseline << ")!" << std::endl;
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

    const int flav = abs(slc->truth.initpdg);
    if (flav != 12 && flav != 14) return;

    const int hcIdx = 0; // TODO TODO always FHC
    const int flavIdx = (flav == 12) ? 0 : 1;
    const int signIdx = (slc->truth.initpdg > 0) ? 0 : 1;

    TH1* h = fScale[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(slc->truth.E);

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
        GetNuMIFluxSyst("fractional_uncertainties/hadron" + name, "hfrac_hadron_", name));
    return ret;
  }

  //----------------------------------------------------------------------
  std::vector<const ISyst*> GetNuMIBeamlineFluxSysts()
  {
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

    std::vector<const ISyst*> ret;
    for (const auto& name : syst_names)
      ret.push_back(GetNuMIFluxSyst("beam_systematic_shifts", "hsyst_beam_", name));
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
    std::vector<const ISyst*> ret = GetNuMIBeamlineFluxSysts();
    std::vector<const ISyst*> add = GetNuMIPCAFluxSysts(Npcs);
    const ISyst* stat = GetNuMIFluxStatisticalUncertainty();
    ret.insert(ret.end(), add.begin(), add.end());
    ret.push_back(stat);
    return ret;
  }
} // namespace ana
