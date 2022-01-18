#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Ratio.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "sbnana/CAFAna/XSec/Flux.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"

using namespace ana;

void test_xsec_flux()
{
  //  SpectrumLoader loader("/icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/ICARUS_prod_2020A_00_numioffaxis_v09_10_01_reco2/ICARUS_prod_2020A_00_numioffaxis_v09_10_01_reco2_CAFMaker_out_flat.root");

  // The only file currently existing with the correct genie inttype filled
  SpectrumLoader loader("/sbnd/app/users/bckhouse/dev/simulation_genie_icarus_numi_volDetEnclosure_20201215T090631_113-0083_gen_20201215T104338_filter_20201215T115238_g4_20201215T215526_detsim_20210113T201733_reco1_20210114T015949_reco2.caf.root");

  const Cut fidvol([](const caf::SRSliceProxy* sr)
                   {
                     const auto& v = sr->truth.position;

                     // Didn't simulate the other cryostat
                     return ((//( v.x < -71.1 - 25 && v.x > -369.33 + 25 ) ||
                              ( v.x > 71.1 + 25 && v.x < 369.33 - 25 )) &&
                             (( v.y > -181.7 + 25 && v.y < 134.8 - 25 ) &&
                              ( v.z > -895.95 + 30 && v.z < 895.95 - 50 ) ));
                   });

  FluxTimesNuclei flux(loader, Binning::Simple(50, 0, 10), fidvol, 14);

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
}
