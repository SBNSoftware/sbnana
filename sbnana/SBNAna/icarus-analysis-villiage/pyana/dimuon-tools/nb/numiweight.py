import uproot
import numpy as np
import pandas as pd

#JD Note June 12, 2025: The ppfx weights are taken care of in numisyst.py script. This script WAS used to fix the concrete issue. Now it will be used to fix concrete AND other CV adjustments needed to reweight from original flux version to latest version. 

## JD: June 12: I commented this out cause I think it's never used (within this script). If nothing is broken when I run it next then all good and can delete.
#FSYST = "/icarus/data/users/gputnam/icarus_numi_flux_syst_ana_v2.root"
#FSYST = "/exp/icarus/data/users/gputnam/thesis_work/icarus_numi_flux_syst_ana_v2.root"
## TODO: Replace above with /exp/icarus/data/users/awood/2024-10-03_out_450.37_7991.98_79512.66.root 

def histdf(h):
    values = h.values()
    bins = h.axis().edges()
    idx = pd.IntervalIndex.from_breaks(bins)
    return pd.Series(values, idx, name=None)
           
def getallpdg_histdf(d, prefix, pdgs=None):
    if pdgs is None:
        pdgs = {
            "numu": 14,
            "numubar": -14,
            "nue": 12,
            "nuebar": -12   
        }

    hs = []
    for pdgname, pdgcode in pdgs.items():
        h = histdf(d[prefix + pdgname])
        h.index = pd.MultiIndex.from_product([[pdgcode], h.index], names=["pdg", "E"])
        hs.append(h)
    return pd.concat(hs)

## JD: June 12: I commented this out cause I think it's never used. If nothing is broken when I run it next then all good and can delete.
#def cv(nupdg, nuE):
#    flux_f = uproot.open(FSYST)
#    wgts = getallpdg_histdf(flux_f["ppfx_flux_weights"], "hweights_fhc_")
#    nuind = pd.MultiIndex.from_arrays([nupdg, nuE])
#    iloc = wgts.index.get_indexer(nuind)
#    match_wgts = wgts.iloc[iloc]
#    match_wgts.loc[iloc < 0, :] = 1.
#    match_wgts.index = nuE.index
#    return match_wgts

# Below is how Gray accounted for concrete winter 2024:
### TODO: get rid of above file and everything below. Was included to change CV according to g3Chase, but now that is just #included in the new file Tony sent me.
#FLUX_CORRECTION_FILE = "/exp/icarus/data/users/awood/numi_flux_weights_g3Chase/g3Chase_weights.root"
#FLUX_CORRECTION_HIST_PREFIX = "hweights_fhc_"
#def concrete_cv(nupdg, nuE):
#    pdgs = {
#        "numu_total": 14,
#        "numubar_total": -14,
#        "nue_total": 12,
#        "nuebar_total": -12, 

#        # use parent PDG to reweight
#        "numu_kpm": 321, 
#        "numubar_kpm": -321, 
#        "nue_k0l": 130
#    }

#    fluxcorr = getallpdg_histdf(uproot.open(FLUX_CORRECTION_FILE), FLUX_CORRECTION_HIST_PREFIX, pdgs)

#    nuind = pd.MultiIndex.from_arrays([nupdg, nuE])
#    iloc = fluxcorr.index.get_indexer(nuind)
#    match_corr = fluxcorr.iloc[iloc]
#    match_corr.loc[iloc < 0, :] = 1.
#    match_corr.index = nupdg.index

#    return match_corr

# stuff below added by Jamie 6/12/25, to replace concrete function above (these reweights include concrete effect and others).

FLUX_CORRECTION_FILE = "/pnfs/icarus/persistent/users/faabdalr/2025-04-08_out_450.37_7991.98_79512.66.root"
FLUX_CORRECTION_FILE = "/exp/icarus/app/users/jdyer/2025-04-08_out_450.37_7991.98_79512.66.root" # This is a copy of above, that I can access.
FLUX_CORRECTION_HIST_PREFIX = "hnom_"

def update_flux_version(nupdg, nuE):
    pdgs = {
        "numu_weights": 14,
        "numubar_weights": -14,
        "nue_weights": 12,
        "nuebar_weights": -12, 
        # or, use parent PDG to reweight:
        "numu_kpm_weights": 321, 
        "numubar_kpm_weights": -321, 
        "nue_k0l_weights": 130
    }
    fluxcorr = getallpdg_histdf(uproot.open(FLUX_CORRECTION_FILE)["g4numi_reweight_v02_00-->v03_02"]["fhc"], FLUX_CORRECTION_HIST_PREFIX, pdgs)
    nuind = pd.MultiIndex.from_arrays([nupdg, nuE])
    iloc = fluxcorr.index.get_indexer(nuind)
    match_corr = fluxcorr.iloc[iloc]
    match_corr.loc[iloc < 0, :] = 1.
    match_corr.index = nupdg.index
    return match_corr
