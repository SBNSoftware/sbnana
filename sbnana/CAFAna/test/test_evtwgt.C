#include "sbnana/CAFAna/Core/Ratio.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SystShifts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
using namespace ana;

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/TruthCuts.h"

#include "TH1.h"
#include "TPad.h"

const double pot = 6e20;

const std::vector<std::string> systs = {
  "DISAth",
  "DISBth",
  "DISCv1u",
  "DISCv2u",
  "IntraNukeNabs",
  "IntraNukeNcex",
  "IntraNukeNinel",
  "IntraNukeNmfp",
  "IntraNukeNpi",
  "IntraNukePIabs",
  "IntraNukePIcex",
  "IntraNukePIinel",
  "IntraNukePImfp",
  "IntraNukePIpi",
  "NC",
  "NonResRvbarp1pi",
  //  "NonResRvbarp1pi", // need to do something with the Alt variants...
  "NonResRvbarp2pi",
  //  "NonResRvbarp2pi",
  "NonResRvp1pi",
  //  "NonResRvp1pi",
  "NonResRvp2pi",
  //  "NonResRvp2pi",
  "ResDecayGamma",
  "CCResAxial",
  "CCResVector",
  //  "CohMA", // return lots of NaNs...
  //  "CohR0",
  "NCELaxial",
  "NCELeta",
  "NCResAxial",
  "NCResVector",
  "QEMA",
};

void test_evtwgt()
{
  SpectrumLoader loader("workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_proton_genie_intrnue_spill_gsimple-configf-v1_tpc_flat_caf_sbnd");

  HistAxis axis("True Energy (GeV)", Binning::Simple(20, 0, 5), SIMPLEVAR(truth.E));

  std::map<std::string, Spectrum*> sPP, sP, sN, sNN;
  for(const std::string& syst: systs){
    ISyst* s = new SBNWeightSyst(syst);
    sPP[syst] = new Spectrum(loader, axis, kIsCC, SystShifts(s, +2));
    sP [syst] = new Spectrum(loader, axis, kIsCC, SystShifts(s, +1));
    sN [syst] = new Spectrum(loader, axis, kIsCC, SystShifts(s, -1));
    sNN[syst] = new Spectrum(loader, axis, kIsCC, SystShifts(s, -2));
  }

  Spectrum snom(loader, axis, kIsCC);

  // This is very slow (to be investigated), so commented out for now
  //  std::vector<Spectrum> multiverse;
  //  multiverse.reserve(100);
  //  for(int i = 0; i < 100; ++i) multiverse.emplace_back(loader, axis, kIsCC, kNoShift, GetUniverseWeight(systs, i));

  loader.Go();

  for(const std::string& syst: systs){
    TH1* h = sPP[syst]->ToTH1(pot, kRed);
    h->SetTitle(syst.c_str());
    h->Draw("hist");
    sP[syst]->ToTH1(pot, kRed-7)->Draw("hist same");
    sN[syst]->ToTH1(pot, kBlue-7)->Draw("hist same");
    sNN[syst]->ToTH1(pot, kBlue)->Draw("hist same");
    snom.ToTH1(pot)->Draw("hist same");
    gPad->Print(("spect_"+syst+".pdf").c_str());

    h = (*sPP[syst] / snom).ToTH1(kRed);
    h->SetTitle(syst.c_str());
    h->GetYaxis()->SetRangeUser(.7, 1.3);
    h->Draw("hist");
    (*sP [syst] / snom).ToTH1(kRed-7)->Draw("hist same");
    (*sN [syst] / snom).ToTH1(kBlue-7)->Draw("hist same");
    (*sNN[syst] / snom).ToTH1(kBlue)->Draw("hist same");
    gPad->Print(("ratio_"+syst+".pdf").c_str());
  }

  // TH1* h = snom.ToTH1(pot);
  // h->Draw("hist");
  // h->GetYaxis()->SetRangeUser(0, 1.3*h->GetMaximum());
  // for(const Spectrum& s: multiverse) s.ToTH1(pot, kGray)->Draw("hist same");
  // h->Draw("hist same");
  // gPad->Print("multiverse.pdf");
}
