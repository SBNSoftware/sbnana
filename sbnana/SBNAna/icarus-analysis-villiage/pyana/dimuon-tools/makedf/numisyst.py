import uproot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#FSYST = "/exp/icarus/data/users/gputnam/thesis-work/icarus_numi_flux_syst_ana_v2.root"
FSYST = "/pnfs/icarus/persistent/users/faabdalr/2025-04-08_out_450.37_7991.98_79512.66.root" # change made June 12, 2025 to update to most recent flux weights from Fatima and Dan. This change takes care of the ppfx weights. (See notes saved to my desktop in muons_analysis folder titled 250611_flux_wgt_update.)
FSYST = "/exp/icarus/app/users/jdyer/2025-04-08_out_450.37_7991.98_79512.66.root" # This is a copy of above, that I can access.

# original
#beam_uncertainties = [
#    "beam_div",
#    "beam_shift_x",
#    "beam_spot",
#    "horn1_x",
#    "horn1_y",
#    "horn_current_plus",
#    "water_layer"
#    # Question 6/20/25: is there a reason Gray didn't include the beam_shift_y ones that are also in the file?
#    # My guessed answer: because all the other ones were applied as symmetric +- variations, but the y ones have specific + and -. 
#] 

# Question: why didn't Gray include hadron production uncertainties? It existed in FSYST file --> fractional_uncertainties --> hadron. Answer: the pca is accounting for the hadron uncertainties and can be used instead of the stuff labeled "hadron" in the file (says Gray.)

# Jamie update 6/20/25:
focusing_uncertainties_names = [ # the +1 sigma variations.
    "Beam_shift_x",
    "Beam_shift_y",
    "Beam_spot",
    "Horn1_x",
    "Horn1_y",
    "Horn2_x",
    "Horn2_y",
    "Horn_current",
    "Horn_water",
    "Target_z"
]

ps1_focusing_uncertainties = [ # the +1 sigma variations.
    "Beam_shift_x_p1mm",
    "Beam_shift_y_p1mm",
    "Beam_spot_1_7mm",
    "Horn1_x_p3mm",
    "Horn1_y_p3mm",
    "Horn2_x_p3mm",
    "Horn2_y_p3mm",
    "Horn_p2kA",
    "Horns_2mm_water",
    "Target_z_p7mm"
]
ms1_focusing_uncertainties = [ # the -1 sigma variations.
    "Beam_shift_x_m1mm",
    "Beam_shift_y_m1mm",
    "Beam_spot_1_3mm",
    "Horn1_x_m3mm",
    "Horn1_y_m3mm",
    "Horn2_x_m3mm",
    "Horn2_y_m3mm",
    "Horn_m2kA",
    "Horns_0mm_water",
    "Target_z_m7mm"
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

def numisyst(nupdg, nuE, fsyst=FSYST, symmetrize_errBars=True):
    flux_f = uproot.open(fsyst)

    cv = getallpdg_histdf(flux_f["ppfx_flux_weights"], "hweights_fhc_")
    cv.name = ("ppfx", "cv")

    beam_syst_wgts = []
    
    ## original:
    #for uc in beam_uncertainties:
    #    uncdf = getallpdg_histdf(flux_f["fractional_uncertainties"]["beam"][uc], "hfrac_beam_" + uc + "_fhc_")
    #    wgtdf_p = 1 + uncdf
    #    wgtdf_m = 1 - uncdf # JD June 30, 2025: I think this is/was the wrong way to do this! But I also don't think wgtdf_m ever got referenced. Nvm, it's fine.
    #    wgtdf_p.name = (uc, "ps1")
    #    wgtdf_m.name = (uc, "ms1")
    #    beam_syst_wgts.append(wgtdf_p)
    #    beam_syst_wgts.append(wgtdf_m)
    
    # Jamie update 6/20/25:
    for iuc, uc in enumerate(focusing_uncertainties_names):
        uncdf_p = getallpdg_histdf(flux_f["beam_focusing_uncertainties"]["fhc"], "hsyst_beam_" + ps1_focusing_uncertainties[iuc] + "_fhc_")
        uncdf_m = getallpdg_histdf(flux_f["beam_focusing_uncertainties"]["fhc"], "hsyst_beam_" + ms1_focusing_uncertainties[iuc] + "_fhc_")
        if symmetrize_errBars:
            uncdf = (np.abs(uncdf_p) + np.abs(uncdf_m))/2. # take the average b/c the covariance matrix requires symmetric uncertainties.
            #(Note: the "plus" and "minus" refer to the knobs, not the direction of event rate change. So distinguising between plus and minus is arbitrary anyway -- average should be totally fine.)
            wgtdf_p = 1 + uncdf
            wgtdf_m = 1 - uncdf # consider the "negative" knob turn to be the opposite sign change to the flux as the "positive" knob turn.
            # also save the actual ones for studying downstream:
            # (The reason I don't do this in same 'if' above is b/c I want the 'ps1' and 'ms1' labels saved first.)
            wgtdf_pActual = 1 + uncdf_p
            wgtdf_mActual = 1 + uncdf_m
            wgtdf_pActual.name = (uc, "actual_ps1")
            wgtdf_mActual.name = (uc, "actual_ms1")
            beam_syst_wgts.append(wgtdf_pActual)
            beam_syst_wgts.append(wgtdf_mActual)
            #print(beam_syst_wgts)
            
        else: # note: if running this else, I need to figure something other than the covariance matrix out.
            wgtdf_p = 1 + uncdf_p
            wgtdf_m = 1 + uncdf_m
        wgtdf_p.name = (uc, "ps1")
        wgtdf_m.name = (uc, "ms1")
        beam_syst_wgts.append(wgtdf_p)
        beam_syst_wgts.append(wgtdf_m)          
    
    for i in pca_components: # note: the PCA components account for the hadron production uncertainties.
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


