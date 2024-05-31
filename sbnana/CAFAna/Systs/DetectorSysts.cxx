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
  int fIdx,
  const std::string &systFile)
    : ISyst("detector_" + name, "Detector Model: " + name)
    , fHistName(name)
    , fName(name)
    , fDualSided(fIdx)
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

      fSystFilePath = std::string(sbndata) + "detectorSysts/detsyst_ratios_combinedsysts.root";
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
    std::string hName = "";
    //enable either 1 histogram or 2 for a +/- 1 sigma type effect, think about improving later
    if(fDualSided == 0)
    {
      hName = "hist_recoNuE_1mNp_syst_" + fHistName;
    }
    //need to make sure that the ratio histograms are named such to allow for this implementation, figure out a smarter way later
    if(fDualSided == 1)
    {
      if(sigma > 0)
      {
        hName = "hist_recoNuE_1mNp_syst_" + fHistName + "_p1sigma";
      }
      if(sigma < 0)
      {
        hName = "hist_recoNuE_1mNp_syst_" + fHistName + "_m1sigma";
      }
    }
    TH1* h = (TH1*)f.Get(hName.c_str());

    double recoNuE = kIcarus202401RecoNuE(slc);
    const int bin = h->FindBin(recoNuE);

    if (bin == 0 || bin == h->GetNbinsX() + 1) return;
    // attempt to not change any weight if the value is inf or nan, just continue with weight as-is... 
    if(fDualSided == 0)
    {                                                    
      std::string s_m1 = "m1sigma";
      bool found_m1 = fHistName.find(s_m1) != std::string::npos; 
      if((std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin))) || h->GetBinContent(bin) == 0)
        weight *= 1;
      else
      {
        if(found_m1)
          weight *= 1 - sigma * (h->GetBinContent(bin)-1.0);
        else
          weight *= 1 + sigma * (h->GetBinContent(bin)-1.0);
      }
    }
    if(fDualSided == 1)
    {
      if ((std::isinf(h->GetBinContent(bin)) || std::isnan(h->GetBinContent(bin))) || h->GetBinContent(bin) == 0)
        weight *= 1;
      else
        weight *= 1 + std::abs(sigma) * (h->GetBinContent(bin)-1.0);
    }
  }   

  const DetectorSyst* GetDetectorSyst(const std::string& dir,
                                      const std::string& prefix,
                                      const std::string& name,
                                      int fIdx,
                                      const std::string& systFile /* = "" */)
  {
    // Make sure we always give the same one back                                                                                                                
    static std::map<std::string, const DetectorSyst*> cache;

    const std::string key = dir + "/" + prefix + name;
    if (cache.count(key) == 0) cache[key] = new DetectorSyst(dir, prefix, name, fIdx, systFile);

    return cache[key];
  }

  std::vector<const ISyst*> GetDetectorSysts()
  {
    //try to keep things in the same order, otherwise it gets very complicated but it should work as is
    
    static std::map<std::string, int> syst_names; 
    /*
    syst_names.insert(std::make_pair("var1_untunedtpcsigshape", 0));
    syst_names.insert(std::make_pair("var2_tpcind2transparency", 1));
    syst_names.insert(std::make_pair("var3_tpcind1gain", 1));
    syst_names.insert(std::make_pair("var4_pmtdecreasedqe", 0));
    syst_names.insert(std::make_pair("var5_ellipsoidalrecomb", 0));
    syst_names.insert(std::make_pair("var6_tpccohnoise", 1));
    syst_names.insert(std::make_pair("var7_tpcintnoise", 1));
    */
    syst_names.insert(std::make_pair("var1_untunedtpcsigshape", 0));
    syst_names.insert(std::make_pair("var2_tpcind2transparency_m1sigma", 0));
    syst_names.insert(std::make_pair("var2_tpcind2transparency_p1sigma", 0));
    syst_names.insert(std::make_pair("var3_tpcind1gain_p1sigma", 0));
    syst_names.insert(std::make_pair("var3_tpcind1gain_m1sigma", 0));
    syst_names.insert(std::make_pair("var4_pmtdecreasedqe", 0));
    syst_names.insert(std::make_pair("var5_ellipsoidalrecomb", 0));
    syst_names.insert(std::make_pair("var6_tpccohnoise_p1sigma", 0));
    syst_names.insert(std::make_pair("var6_tpccohnoise_m1sigma", 0));
    syst_names.insert(std::make_pair("var7_tpcintnoise_p1sigma", 0));
    syst_names.insert(std::make_pair("var7_tpcintnoise_m1sigma", 0));
    syst_names.insert(std::make_pair("var8_nullvar", 0));
    static std::vector<const ISyst*> ret;
    if(!ret.empty()) return ret;

    for(auto [name, mode]: syst_names)
      ret.push_back(GetDetectorSyst("detector_systematic_shifts", "hist_recoNuE_1mNp_syst_" , name, mode));

    return ret;
  }
}
