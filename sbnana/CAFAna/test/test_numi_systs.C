#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnana/CAFAna/Systs/NuMIFluxSysts.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"

using namespace ana;

void test_numi_systs()
{
  const std::string fname = "/icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/ICARUS_prod_2020A_00_numioffaxis_v09_10_01_reco2/ICARUS_prod_2020A_00_numioffaxis_v09_10_01_reco2_CAFMaker_out_flat.root";

  const std::vector<const ISyst*> systs = GetAllNuMIFluxSysts(5);

  SpectrumLoader loader(fname);

  const Var kTrueE = SIMPLEVAR(truth.E);

  HistAxis axis("True neutrino energy (GeV)",
                Binning::Simple(50, 0, 5),
                kTrueE);

  Spectrum snom(loader, axis, kNoCut);
  std::vector<Spectrum*> sup, sdn;

  for(const ISyst* s: systs){
    sup.push_back(new Spectrum(loader, axis, kNoCut, SystShifts(s, +1)));
    sdn.push_back(new Spectrum(loader, axis, kNoCut, SystShifts(s, -1)));
  }

  loader.Go();

  const double pot = 1e20;
  TH1* hnom = snom.ToTH1(pot);

  for(unsigned int i = 0; i < systs.size(); ++i){
    new TCanvas;
    TH1* h = (TH1*)hnom->Clone(UniqueName().c_str());
    h->SetTitle(systs[i]->LatexName().c_str());
    h->Draw("hist");
    sup[i]->ToTH1(pot, kRed)->Draw("hist same");
    sdn[i]->ToTH1(pot, kBlue)->Draw("hist same");
    gPad->Print((systs[i]->ShortName()+".png").c_str());
  }
}
