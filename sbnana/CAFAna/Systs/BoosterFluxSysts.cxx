#include "sbnana/CAFAna/Systs/BoosterFluxSysts.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>

namespace ana {

//----------------------------------------------------------------------
BoosterFluxHadronSyst::~BoosterFluxHadronSyst()
{
  for(int i = 0; i < 3; ++i) delete fScale[i];
}

//----------------------------------------------------------------------
void BoosterFluxHadronSyst::Shift(double sigma, caf::SRSliceProxy* slc, double &weight) const
{
  // Only implemented for numus so far
  if(/*sr->mc.nnu == 0 ||*/ slc->truth.initpdg != 14) return;

  if(!fScale[0]){
    TFile f((FindCAFAnaDir() + "/Systs/booster_flux_shifts.root").c_str());
    assert(!f.IsZombie());

    for(int det = 0; det < 3; ++det){
      std::string detStr;
      if(det == 0) detStr = "nd";
      if(det == 1) detStr = "ub";
      if(det == 2) detStr = "fd";

      TH1*& h = fScale[det];

      h = (TH1*)f.Get(TString::Format("syst%d/%s_numu", fIdx, detStr.c_str()).Data());
      h = (TH1*)h->Clone(UniqueName().c_str());
      assert(h);
      h->SetDirectory(0);
    }
  } // end if

  // TODO use the regular detector flag as soon as it's filled reliably
  int det;
  const double L = slc->truth.baseline;
  if(L < 150) det = 0;
  else if(L < 500) det = 1;
  else det = 2;

  TH1* h = fScale[det];
  assert(h);

  const int bin = h->FindBin(slc->truth.E);

  if(bin == 0 || bin == h->GetNbinsX()+1) return;

  weight *= exp(sigma*h->GetBinContent(bin));
}

//----------------------------------------------------------------------
const BoosterFluxHadronSyst* GetBoosterFluxHadronSyst(unsigned int i)
{
  // Make sure we always give the same one back
  static std::vector<const BoosterFluxHadronSyst*> cache;

  if(i >= cache.size()) cache.resize(i+1);

  if(!cache[i]) cache[i] = new BoosterFluxHadronSyst(i);

  return cache[i];
}

//----------------------------------------------------------------------
BoosterFluxHadronSystVector GetBoosterFluxHadronSysts(unsigned int N)
{
  BoosterFluxHadronSystVector ret;
  for(unsigned int i = 0; i < N; ++i)
    ret.push_back(GetBoosterFluxHadronSyst(i));

  return ret;
}

} // namespace ana
