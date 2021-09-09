#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>

namespace ana
{

  //----------------------------------------------------------------------
  NuMIFluxSyst::~NuMIFluxSyst()
  {
    for(int i = 0; i < 2; ++i)
      for(int j = 0; j < 2; ++j)
        for(int k = 0; k < 2; ++k)
          delete fScale[i][j][k];
  }

  //----------------------------------------------------------------------
  void NuMIFluxSyst::Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const
  {
    if(slc->truth.index < 0) return; // No neutrino

    // TODO switch to regular detector flag as soon as it's filled reliably
    if(slc->truth.baseline < 500){
      std::cout << "Using NuMIFluxSyst on what appears not to be an icarus file (baseline = " << slc->truth.baseline << ")!" << std::endl;
      abort();
    }

    if(!fScale[0][0][0]){
      // TODO - where will this file live / how will it be found?
      TFile f("icarus_numi_flux_syst_ana.root");
      assert(!f.IsZombie());

      for(int hcIdx: {0, 1}){
        for(int flavIdx: {0, 1}){
          for(int signIdx: {0, 1}){
            std::string hName = "fractional_uncertainties/hfrac_"+fName;
            if(hcIdx == 0) hName += "_fhc"; else hName += "_rhc";
            if(flavIdx == 0) hName += "_nue"; else hName += "_numu";
            if(signIdx == 1) hName += "bar";

            TH1* h = (TH1*)f.Get(hName.c_str());
            assert(h);
            h = (TH1*)h->Clone(UniqueName().c_str());
            h->SetDirectory(0);

            fScale[hcIdx][flavIdx][signIdx] = h;
          }
        }
      }
    } // end if

    const int flav = abs(slc->truth.initpdg);
    if(flav != 12 && flav != 14) return;

    const int hcIdx = 0; // TODO TODO always FHC
    const int flavIdx = (flav == 12) ? 0 : 1;
    const int signIdx = (slc->truth.initpdg > 0) ? 0 : 1;

    TH1* h = fScale[hcIdx][flavIdx][signIdx];
    assert(h);

    const int bin = h->FindBin(slc->truth.E);

    if(bin == 0 || bin == h->GetNbinsX()+1) return;

    weight *= 1 + sigma*h->GetBinContent(bin);
  }

  //----------------------------------------------------------------------
  const NuMIFluxSyst* GetNuMIFluxSyst(const std::string& name)
  {
    // Make sure we always give the same one back
    static std::map<std::string, const NuMIFluxSyst*> cache;

    if(cache.count(name) == 0) cache[name] = new NuMIFluxSyst(name);

    return cache[name];
  }

  //----------------------------------------------------------------------
  std::vector<const ISyst*> GetNuMIFluxSysts()
  {
    const std::vector<std::string> syst_names = {"thintarget", "pCpi", "pCk", "pCnu", "nCpi", "mesinc", "nua", "nuAlFe", "attenuation", "others"};

    std::vector<const ISyst*> ret;
    for(std::string name: syst_names) ret.push_back(GetNuMIFluxSyst(name));
    return ret;
  }

} // namespace ana
