#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "cafanacore/Ratio.h"
#include "cafanacore/Spectrum.h"

#include "sbnana/CAFAna/XSec/Flux.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TPad.h"

#include <iostream>

using namespace ana;

void test_xsec_flux()
{
  // Need to export SAM_EXPERIMENT=sbn
  SpectrumLoader loader("Official_icarusprod_2021C_NUMI_Nu_Cosmics_v09_37_01_03p01_caf");

  const NuTruthCut fidvol([](const caf::SRTrueInteractionProxy* nu)
                   {
                     const auto& v = nu->position;

                     // Didn't simulate the other cryostat
                     return ((//( v.x < -71.1 - 25 && v.x > -369.33 + 25 ) ||
                              ( v.x > 71.1 + 25 && v.x < 369.33 - 25 )) &&
                             (( v.y > -181.7 + 25 && v.y < 134.8 - 25 ) &&
                              ( v.z > -895.95 + 30 && v.z < 895.95 - 50 ) ));
                   });

  FluxTimesNuclei flux(loader.NuTruths(), Binning::Simple(50, 0, 5), fidvol, 14);

  loader.Go();

  const double V = (((369.33 - 25) - (71.1 + 25)) *
                    ((134.8 - 25 ) - ( -181.7 + 25)) *
                    ((895.95 - 50 ) -( -895.95 + 30))) * 1e-6; // cubic metres

  const double m = V*1.3982*1000; // kg

  const double amu = 1.66e-27;

  const double N = m / (39.948*amu); // number of nuclei

  std::cout << "V = " << V << ", m = " << m << " N = " << N << std::endl;

  new TCanvas;
  flux.ToTH1(1/N, kBlack, kSolid, kBinDensity)->Draw("hist");

  //  new TCanvas;
  TFile* fin = new TFile("/sbnd/app/users/bckhouse/dev/numi-at-icarus-flux-systematics/icarus_numi_flux_syst_ana.root");
  if(!fin->IsZombie()){
    TH1* hnom = (TH1*)fin->Get("ppfx_output/run15/fhc/nom/hnom_numu");
    hnom->SetLineColor(kRed);
    hnom->Draw("hist same");
  }

  gPad->Print("test_xsec_flux.pdf");
}
