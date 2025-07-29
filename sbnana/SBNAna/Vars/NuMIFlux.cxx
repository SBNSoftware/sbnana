#include "sbnana/SBNAna/Vars/NuMIFlux.h"

#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ana {

  //===============================================================
  // 05/19/25 JK) This one uses "2025-04-08_out_450.37_7991.98_79512.66.root",
  //              to provide a CV correction of beam width setting 1.5mm/1.4mm;
  //              NuMI2023 reprocessing simulation was done with CV beam width of 1.4mm,
  //              but the data is more like 1.5mm.
  NuMIFluxCorrection::NuMIFluxCorrection()
  { 
      
    std::cout << "[NuMIFluxCorrection::NuMIFluxCorrection] Called" << std::endl;

    const char* sbndata = std::getenv("SBNDATA_DIR");
    if (!sbndata) {
      std::cout << "NuMIFluxCorrection: $SBNDATA_DIR environment variable not set. Please setup "
                   "the sbndata product."
                << std::endl;
      std::abort(); 
    } 
        
    // normal PPFX weights
    fFluxFilePath = std::string(sbndata) +
                   "beamData/NuMIdata/2025-04-08_out_450.37_7991.98_79512.66.root";

    TFile f(fFluxFilePath.c_str());
    if (f.IsZombie()) {
      std::cout << "NuMIFluxCorrection: Failed to open " << fFluxFilePath << std::endl;
      std::abort();
    }

    for (int hcIdx : {0, 1}) {
      for (int flavIdx : {0, 1}) {
        for (int signIdx : {0, 1}) {
          std::string hNamePPFX = "ppfx_flux_weights/hweights_";
          if (hcIdx == 0)
            hNamePPFX += "fhc_";
          else
            hNamePPFX += "rhc_";
          if (flavIdx == 0)
            hNamePPFX += "nue";
          else
            hNamePPFX += "numu";
          if (signIdx == 1) hNamePPFX += "bar";
    
          TH1* h_ppfx = (TH1*)f.Get(hNamePPFX.c_str());
          if (!h_ppfx) {
            std::cout << "[NuMIFluxCorrection::NuMIFluxCorrection] Failed to find " << hNamePPFX << " from " << f.GetName()
                      << std::endl;
            std::abort();
          }
          h_ppfx = (TH1*)h_ppfx->Clone(UniqueName().c_str());
          h_ppfx->SetDirectory(0);

          fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;
        }
      }
    }

    // G3Chase+G4Update+BeamWidth(1.4mm->1.5mm)

    for (int currIdx : {0, 1}) {
      for (int flavIdx : {0, 1}) {
        for (int signIdx : {0, 1}) {
          for (int pdgIdx : {0, 1, 2, 3}) {

            std::string hNameCVCorr = "g4numi_reweight_v02_00-->v03_02";

            // hnom_nue_k0l_weights

            // horn current
            if(currIdx==0) hNameCVCorr += "/fhc";
            else           hNameCVCorr += "/rhc";

            // nu flavor
            if(flavIdx == 0) hNameCVCorr += "/hnom_nue";
            else             hNameCVCorr += "/hnom_numu";

            // nu/anti-nu
            if(signIdx == 1) hNameCVCorr += "bar";

            // parent particle type
            if(pdgIdx==0)      hNameCVCorr += "_pipm";
            else if(pdgIdx==1) hNameCVCorr += "_kpm";
            else if(pdgIdx==2) hNameCVCorr += "_mu";
            else               hNameCVCorr += "_k0l";

            hNameCVCorr += "_weights";

            TH1* h_CVCorr = (TH1*)f.Get(hNameCVCorr.c_str());
            if (!h_CVCorr) {
              // We don't expect nue(bar)_pipm
              if( flavIdx==0 && pdgIdx==0 ){
                printf("[NuMIFluxCorrection::NuMIFluxCorrection] Skip: (currIdx, flavIdx, signIdx, pdgIdx) = (%d, %d, %d, %d)\n", currIdx, flavIdx, signIdx, pdgIdx);
                fWeightCVCorr[currIdx][flavIdx][signIdx][pdgIdx] = nullptr;
                continue;
              }
              printf("[NuMIFluxCorrection::NuMIFluxCorrection] Fail: (currIdx, flavIdx, signIdx, pdgIdx) = (%d, %d, %d, %d)\n", currIdx, flavIdx, signIdx, pdgIdx);            
              printf("[NuMIFluxCorrection::NuMIFluxCorrection]   %s\n", hNameCVCorr.c_str());
              abort();
            }
            else{
              printf("[NuMIFluxCorrection::NuMIFluxCorrection] Found: (currIdx, flavIdx, signIdx, pdgIdx) = (%d, %d, %d, %d)\n", currIdx, flavIdx, signIdx, pdgIdx);
              printf("[NuMIFluxCorrection::NuMIFluxCorrection]   %s\n", hNameCVCorr.c_str());
            }
            h_CVCorr = (TH1*)h_CVCorr->Clone(UniqueName().c_str());
            h_CVCorr->SetDirectory(0);

            fWeightCVCorr[currIdx][flavIdx][signIdx][pdgIdx] = h_CVCorr;
          } // END loop PDG index
        } // END loop nu type
      } // END loop nu flavor
    } //END loop horn current
  }

  NuMIFluxCorrection& NuMIFluxCorrection::Instance()
  {
    static NuMIFluxCorrection m;
    return m;
  }

  double NuMIFluxCorrection::GetWeightFromSRTrueInt(const caf::SRTrueInteractionProxy* nu) const
  {
    if (nu->index < 0 || abs(nu->initpdg) == 16) return 1.0;

    if (!fWeight[0][0][0] || !fWeightCVCorr[0][0][0][1]) {
      std::cout << "Trying to access un-available weight array..." << std::endl;
      std::abort();
    }

    unsigned int hcIdx = 0; // assume always FHC for now...
    unsigned int flavIdx = (abs(nu->initpdg) == 12) ? 0 : 1;
    unsigned int signIdx = (nu->initpdg > 0) ? 0 : 1;
    unsigned int pdgIdx = ParentPDGToIdx(nu->parent_pdg);

    TH1* h = fWeight[hcIdx][flavIdx][signIdx];
    assert(h);

    double weight = 1.0;

    const int bin = h->FindBin(nu->E);
    if ( bin != 0 && bin != h->GetNbinsX() + 1 && !std::isinf(h->GetBinContent(bin)) && !std::isnan(h->GetBinContent(bin)) )
      weight*=h->GetBinContent(bin);

    // parent not applicable
    if ( pdgIdx == 4 ) return weight;

    TH1* h2 = fWeightCVCorr[hcIdx][flavIdx][signIdx][pdgIdx];
    if(!h2){
      return weight;
    }

    double this_NuE = nu->E;
    // muon parent seems to have low stat above 2.0 GeV
    if( pdgIdx==2 ){
      this_NuE = 1.99;
    }

    const int bin2 = h2->FindBin(this_NuE);
    if ( bin2 != 0 && bin2 != h2->GetNbinsX() + 1 && !std::isinf(h2->GetBinContent(bin2)) && !std::isnan(h2->GetBinContent(bin2)) ) {
      weight*=h2->GetBinContent(bin2);
    }

    return weight;
  }

  NuMIFluxCorrection::~NuMIFluxCorrection()
  {
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          delete fWeight[i][j][k];

          for (int l = 0; l < 4; ++k) {

            delete fWeightCVCorr[i][j][k][l];

          }

        }
      }
    }
  }

  unsigned int NuMIFluxCorrection::ParentPDGToIdx(int pdg) const
  {
    if      ( abs(pdg) == 211 ) return 0;
    else if ( abs(pdg) == 321 ) return 1;
    else if ( abs(pdg) == 13  ) return 2;
    else if ( abs(pdg) == 130 ) return 3;
    return 4;
  }

  const TruthVar kGetTruthNuMIFluxCorrection([](const caf::SRTrueInteractionProxy* nu) -> double {
    const NuMIFluxCorrection& m = NuMIFluxCorrection::Instance();
    return m.GetWeightFromSRTrueInt(nu);
  });
  const Var kGetNuMIFluxCorrection([](const caf::SRSliceProxy* slc) -> double {
    return kGetTruthNuMIFluxCorrection(&slc->truth);
  });



}

