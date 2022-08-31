#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/OscillatableSpectrum.h"
#include "cafanacore/Ratio.h"
#include "sbnana/CAFAna/Core/OscCalcSterileApprox.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"

const Var kTrueE = SIMPLEVAR(truth.E);
const Var kRecoE = SIMPLEVAR(reco.reco_energy);

void osc_prob()
{
  const std::string dir = "/sbnd/data/users/jlarkin/workshop_samples/";
  const std::string fnameBeam = dir + "output_SBNOsc_NumuSelection_Modern_SBND.flat.root";

  SpectrumLoader loader(fnameBeam);

  OscillatableSpectrum s(loader.Slices()[kIsNumuCC], HistAxis("True E (GeV)", Binning::Simple(100, .5, 1.5), kTrueE));
  OscillatableSpectrum s_reco(loader.Slices()[kIsNumuCC], HistAxis("Reco E (GeV)", Binning::Simple(100, .5, 1.5), kRecoE));

  loader.Go();

  OscCalcSterileApprox calc;
  calc.SetDmsq(50);
  calc.SetSinSq2ThetaMuMu(1e-2);
  calc.SetSinSq2ThetaMuE(0);

  s.ToTH2(1e20)->Draw("col");
  gPad->Print("osc_prob_2d.png");
  new TCanvas;
  s_reco.ToTH2(1e20)->Draw("col");
  gPad->Print("osc_prob_2d_reco.png");

  Ratio r(s.Oscillated(&calc, 14, 14), s.Unoscillated());
  Ratio r_reco(s_reco.Oscillated(&calc, 14, 14), s_reco.Unoscillated());

  new TCanvas;
  TH1* h = r.ToTH1();
  h->GetYaxis()->SetRangeUser(.99, 1);
  h->GetYaxis()->SetTitle("Ratio to unoscillated");
  h->Draw();

  gPad->Print("osc_prob.png");
  gPad->Print("osc_prob.pdf");
  TFile* fout = new TFile("osc_prob.root", "RECREATE");
  h->Write("osc_prob");

  new TCanvas;
  h = r_reco.ToTH1();
  h->GetYaxis()->SetRangeUser(.99, 1);
  h->GetYaxis()->SetTitle("Ratio to unoscillated");
  h->Draw();

  gPad->Print("osc_prob_reco.png");
  gPad->Print("osc_prob_reco.pdf");

  h->Write("osc_prob_reco");
}
