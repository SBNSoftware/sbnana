#include "sbnana/CAFAna/Systs/DetectorSysts.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Vars/NumuVarsIcarus202401.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ana 
{
  DetectorSyst::DetectorSyst(const std::string &dir,
  const std::string &prefix,
  const std::string &name,
  const std::string &systFile)
    : ISyst("detector_" + name, "Detector Model: " + name)
    , fHistName(name)
    , fName(name)
  {
    if (!systFile.empty()) { fSystFilePath = systFile; }
    else {
      const char* sbndata = std::getenv("SBNDATA_DIR");
      if (!sbndata) {
        std::cout << "DetectorModelSystematics: $SBNDATA_DIR environment variable not set. Please setup "
                     "the sbndata product."
                  << std::endl;
        std::abort();
      }

      fSystFilePath = std::string(sbndata) + "detectorSysts/detector_syst_ratios.root";
      //fSystFilePath = "/exp/icarus/app/users/jzettle/Systematics_SelAnl/Selection/initial_syst_ratios.root";
    }
  }
  void DetectorSyst::Shift(double sigma, caf::SRSliceProxy* slc, double &weight) const
  {
    if (slc->truth.index < 0) return; //No neutrino

    TFile f(fSystFilePath.c_str());
    if (f.IsZombie()) {
        std::cout << "DetectorModelSysts: Failed to open " << fSystFilePath << std::endl;
        abort();
    }
    std::string hName = "ratio_cv_syst_" + fHistName;
    TH1* h = (TH1*)f.Get(hName.c_str());

    double recoNuE = kIcarus202401RecoNuE(slc);
    const int bin = h->FindBin(recoNuE);

    if (bin == 0 || bin == h->GetNbinsX() + 1) return;
    // attempt to not change any weight if the value is inf or nan, just continue with weight as-is...                                                     
    if (std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin)))
      weight *= 1;
    else
      weight *= 1 + sigma * (h->GetBinContent(bin)-1.0);
    
  }   

  const DetectorSyst* GetDetectorSyst(const std::string& dir,
                                      const std::string& prefix,
                                      const std::string& name,
                                      const std::string& systFile /* = "" */)
  {
    // Make sure we always give the same one back                                                                                                                
    static std::map<std::string, const DetectorSyst*> cache;

    const std::string key = dir + "/" + prefix + name;
    if (cache.count(key) == 0) cache[key] = new DetectorSyst(dir, prefix, name, systFile);

    return cache[key];
  }

  std::vector<const ISyst*> GetDetectorSysts()
  {
    const std::vector<std::string> syst_names = {"cohnoise_p1sigma",
                                                 "cohnoise_m1sigma",
                                                 "cohnoise_p2sigma",
                                                 "cohnoise_m2sigma"};
    static std::vector<const ISyst*> ret;
    if(!ret.empty()) return ret;

    for(const std::string& name: syst_names)
      ret.push_back(GetDetectorSyst("detector_systematic_shifts", "hsyst_detector_" , name));

    return ret;
  }
}