import uproot
import numpy as np
import pandas as pd

FSYST = "/icarus/data/users/gputnam/icarus_numi_flux_syst_ana_v2.root"

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

def cv(nupdg, nuE):
    flux_f = uproot.open(FSYST)
    wgts = getallpdg_histdf(flux_f["ppfx_flux_weights"], "hweights_fhc_")
    nuind = pd.MultiIndex.from_arrays([nupdg, nuE])
    iloc = wgts.index.get_indexer(nuind)
    match_wgts = wgts.iloc[iloc]
    match_wgts.loc[iloc < 0, :] = 1.
    match_wgts.index = nuE.index

    return match_wgts

FLUX_CORRECTION_FILE = "/exp/icarus/data/users/awood/numi_flux_weights_g3Chase/g3Chase_weights.root"
FLUX_CORRECTION_HIST_PREFIX = "hweights_fhc_"
def concrete_cv(nupdg, nuE):
    pdgs = {
        "numu_total": 14,
        "numubar_total": -14,
        "nue_total": 12,
        "nuebar_total": -12, 

        # use parent PDG to reweight
	"numu_kpm": 321, 
	"numubar_kpm": -321, 
	"nue_k0l": 130
    }

    fluxcorr = getallpdg_histdf(uproot.open(FLUX_CORRECTION_FILE), FLUX_CORRECTION_HIST_PREFIX, pdgs)

    nuind = pd.MultiIndex.from_arrays([nupdg, nuE])
    iloc = fluxcorr.index.get_indexer(nuind)
    match_corr = fluxcorr.iloc[iloc]
    match_corr.loc[iloc < 0, :] = 1.
    match_corr.index = nupdg.index

    return match_corr
