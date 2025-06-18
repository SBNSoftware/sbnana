import uproot
import numpy as np
import pandas as pd

#FSYST = "/exp/icarus/data/users/gputnam/thesis-work/icarus_numi_flux_syst_ana_v2.root"
FSYST = "/pnfs/icarus/persistent/users/faabdalr/2025-04-08_out_450.37_7991.98_79512.66.root" # change made June 12, 2025 to update to most recent flux weights from Fatima and Dan. This change takes care of the ppfx weights. (See notes saved to my desktop in muons_analysis folder titled 250611_flux_wgt_update.)
FSYST = "/exp/icarus/app/users/jdyer/2025-04-08_out_450.37_7991.98_79512.66.root" # This is a copy of above, that I can access.

beam_uncertainties = [
    "beam_div",
    "beam_shift_x",
    "beam_spot",
    "horn1_x",
    "horn1_y",
    "horn_current_plus",
    "water_layer"
]

pca_components = list(range(20))

def histdf(h):
    values = h.values()
    bins = h.axis().edges()
    idx = pd.IntervalIndex.from_breaks(bins)
    return pd.Series(values, idx, name=None)
           
def getallpdg_histdf(d, prefix):
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

def numisyst(nupdg, nuE, fsyst=FSYST):
    flux_f = uproot.open(fsyst)

    cv = getallpdg_histdf(flux_f["ppfx_flux_weights"], "hweights_fhc_")
    cv.name = ("ppfx", "cv")

    beam_syst_wgts = []
    for uc in beam_uncertainties:
        uncdf = getallpdg_histdf(flux_f["fractional_uncertainties"]["beam"][uc], "hfrac_beam_" + uc + "_fhc_")
        wgtdf_p = 1 + uncdf
        wgtdf_m = 1 - uncdf
        wgtdf_p.name = (uc, "ps1")
        wgtdf_m.name = (uc, "ms1")
        beam_syst_wgts.append(wgtdf_p)
        beam_syst_wgts.append(wgtdf_m)
    
    for i in pca_components:
        uncdf = getallpdg_histdf(flux_f["pca"]["principal_components"], "hpc_%i_fhc_" % i)
        wgtdf_p = 1 + uncdf
        wgtdf_m = 1 - uncdf
        wgtdf_p.name = (("pca%i" % i), "ps1")
        wgtdf_m.name = (("pca%i" % i), "ms1")
        beam_syst_wgts.append(wgtdf_p)
        beam_syst_wgts.append(wgtdf_m)

    wgts = pd.DataFrame([cv] + beam_syst_wgts).T
    nuind = pd.MultiIndex.from_arrays([nupdg, nuE])
    iloc = wgts.index.get_indexer(nuind)
    match_wgts = wgts.iloc[iloc]
    match_wgts.loc[iloc < 0, :] = 1.
    match_wgts.index = nupdg.index 

    return match_wgts


