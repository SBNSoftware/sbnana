import numpy as np
import pandas as pd
import math
import uproot
import random

from pyanalib.dataset import Dataset
from pyanalib.panda_helpers import *
import weights
import reweight_coh
import numiweight
from alp_flux_weights import *
import track_splitting

def simple_dataset(f, key):
    df = pd.read_hdf(f, key=key)
    cv = pd.Series(1, index=df.index)
    cv.name = ("wgt", "cv")
    df = multicol_add(df, cv)
    return Dataset(df, -1, -1)

def simple_dataset_df(df):
    cv = pd.Series(1, index=df.index)
    cv.name = ("wgt", "cv")
    df = multicol_add(df, cv)
    return Dataset(df, -1, -1)

def add_weights(df, wgtdf, tmatch_col=("slc", "tmatch", "idx")):
    na_set = {}
    for w in wgtdf.columns:
        na_set[w] = 1 # don't weight non-truth matched events

    return multicol_merge(df, wgtdf, how="left", left_on=wgtdf.index.names[:2] + [tmatch_col], right_index=True).fillna(na_set)

def mc_dataset(f, key, hdrkey="hdr", mcnukey="mcnuwgt", syst_weights=True, mccut=None, mccut_any=True, isscalar=False, alp=False):
    hdrdf = pd.read_hdf(f, key=hdrkey)
    livetime = 0.
    pot = hdrdf.pot.sum()
    df = pd.read_hdf(f, key)

    mcdf = pd.read_hdf(f, key=mcnukey)
    # fix number of levels
    if mcdf.index.nlevels == 2:
        mcdf["inu"] = mcdf.groupby(level=[0,1]).cumcount()
        mcdf.set_index("inu", append=True, inplace=True)

    # apply any cut on the mc df
    if mccut:
        mcdfcut = mccut(mcdf).groupby(level=[0,1]).any() if mccut_any else mccut(mcdf).groupby(level=[0,1]).all()

        df["mccut"] = mcdfcut
        df = df[df.mccut]

    # Central-Value (just the ppfx), if we can
    try: 
        cv = math.prod([mcdf[w] for w in weights.cv])
        cv.name = ("wgt", "cv", "tot") # just ppfx so far, but will get multiplied by the other cv reweight below.
        df = add_weights(df, pd.DataFrame(cv))
        df[("wgt", "cv", "ppfx", "", "", "")] = df[("wgt", "cv", "tot", "", "", "")]
    except:
        if syst_weights:
            raise
        else:
            df = multicol_add(df, pd.Series(1, index=df.index, name=("wgt", "cv", "tot")))

    # load other corrections (i.e. concrete, correction for horn 1 and 2 mother volume overlap, and reweighting to a 1.5mm beam spot size -- collectively call this "cv other" (as in not ppfx.)) 
    nuE = mcdf.E

    # apply the weight differently for nuetrinos and scalars
    # Neutrinos: use the "total" for that pdg
    # Scalars: lookup the weight for the corresponding parent pdg
    # re-weight the 
    nupdg = mcdf.parent_pdg if isscalar else mcdf.pdg

    #fluxcorr_wgt = numiweight.concrete_cv(nupdg, nuE) # JD 6/12/25
    #fluxcorr_wgt.name = ("wgt", "concrete", "cv", "", "", "") # JD 6/12/25
    #df = add_weights(df, pd.DataFrame(fluxcorr_wgt)) # JD 6/12/25
    #df[("wgt", "cv", "", "", "", "")] *= df[("wgt", "concrete", "cv", "", "", "")] # JD 6/12/25
    fluxcorr_wgt = numiweight.update_flux_version(nupdg, nuE) # JD 6/12/25
    fluxcorr_wgt.name = ("wgt", "cv", "flux_other_than_ppfx", "", "", "") # JD 6/12/25
    df = add_weights(df, pd.DataFrame(fluxcorr_wgt)) # JD 6/12/25
    # JB 04/27/26: Adding calculation of track splitting reweighting here
    cathode_cross_weight, cathode_cross_weight_err, gap_cross_weight, gap_cross_weight_err = track_splitting.correct(df)
    cathode_cross_weight.name = ("wgt", "cv", "cathode_crossing", "", "", "")
    gap_cross_weight.name = ("wgt", "cv", "z_crossing", "", "", "")
    df = add_weights(df, pd.DataFrame(cathode_cross_weight))
    df = add_weights(df, pd.DataFrame(gap_cross_weight))
    df[("wgt", "cv", "tot", "", "", "")] *= df[("wgt", "cv", "flux_other_than_ppfx", "", "", "")]
    df[("wgt", "cv", "tot", "", "", "")] *= df[("wgt", "cv", "cathode_crossing", "", "", "")]
    df[("wgt", "cv", "tot", "", "", "")] *= df[("wgt", "cv", "z_crossing", "", "", "")]
    
    # 12/17/25: Include a reweight for consideration of secondaries for ALP flux
    if alp:
        mch_df = pd.read_hdf(f, key="mch")
        bsm_E = np.array(pd.merge(df, mch_df.E, on=['__ntuple', 'entry'], how='left').E)
        bsm_M = np.array(pd.merge(df, mch_df.M, on=['__ntuple', 'entry'], how='left').M)
        alp_rw = [ALP_flux_rw_for_secondaries(alp_E, (m)) for alp_E, m in zip(bsm_E, bsm_M)]
    else:
        alp_rw = [1]*df.shape[0]
    df[("wgt", "cv", "flux_alp_secondaries", "", "", "")] = alp_rw
    df[("wgt", "cv", "tot", "", "", "")] *= df[("wgt", "cv", "flux_alp_secondaries", "", "", "")]


    print("Generating Coh-weights!")
    # generate coherent pion weights
    cohweight = reweight_coh.cohweight(mcdf, douniv=syst_weights)
    cohweight.columns = pd.MultiIndex.from_tuples([("wgt", "cv", "coh", "", "", "")] + 
                                                  [("wgt", "coh", "univ_%i" % i, "", "", "") for i in range(len(cohweight.columns) - 1)])
    df = add_weights(df, cohweight)
    df[("wgt", "cv", "tot", "", "", "")] *= df[("wgt", "cv", "coh", "", "", "")]
    print("Generated")

    if not syst_weights:
        return Dataset(df, livetime, pot, hdrdf)

    np.random.seed(24601) # My name is Jean-Valjean

    # Systematics
    NUNI = 100
    def make_unidf(systematics, colprefix):
        uni_columns = ["univ_%i" % i for i in range(NUNI)]
        unidf = pd.DataFrame(1, index=mcdf.index, columns=uni_columns)
        for s in systematics:
            print(s)
            if s not in mcdf.columns: 
                print('Cant find systematic called %a in df. Fix this!!' %s)
                continue

            # +/-1 sigma
            if 'ps' in mcdf[s].columns[0][0]: #.startswith("ps"): # JD October 30, 2025
                #ps = mcdf[s].columns[0][0]
                ## ps = "ps1"
                #shift = np.random.normal(size=NUNI) # mean=0, sd=1
                #w = 1 + pd.DataFrame(np.outer((mcdf[s][ps]-1), shift), index=mcdf.index, columns=uni_columns)
                ## JD Todo: average ps and ms using the absolute value difference from one. Done: I did this in numisyst. Okay as long as not super lop-sided.
                ## This seems reasonable if the knob turns affect the event rate in opposite directions.
                ## check how often this is true.
                ## when this is not true, consider treating those with the "morph" treatment below, or some other way of handling the asymmetry.
                
                # November 13, 2025 TODO: 
                #   Change the method. Instead: For +- 1 sigma, randomly choose whether to apply the ps1 or ms1 effect, and apply the chosen one as a one-sided shift. The sign of the chosen one will affect the CV in the correct direction, and the +- effect should average out correctly over the many events. 
                ms1_or_ps1 = random.randint(0,1)
                ps = mcdf[s].columns[ms1_or_ps1][0]
                shift = np.random.normal(size=NUNI)
                w = 1 + pd.DataFrame(np.outer((mcdf[s][ps]-1), np.abs(shift)), index=mcdf.index, columns=uni_columns)
            
        
            # One-sided
            elif mcdf[s].columns[0][0] == "morph":
                shift = np.random.normal(size=NUNI)
                #w = 1 + pd.DataFrame(np.outer((mcdf[s].morph-1)*2, np.abs(shift)), index=mcdf.index, columns=uni_columns)  # I got rid the of the factor of 2 on November 13, 2025. Why should it be there?? I don't think it should be!
                w = 1 + pd.DataFrame(np.outer((mcdf[s].morph-1), np.abs(shift)), index=mcdf.index, columns=uni_columns) 
        
            # Universe uncertainties
            elif mcdf[s].columns[0][0] == "univ_0":
                w = mcdf[s][uni_columns]
                w.columns = [c[0] for c in w.columns]

            unidf = unidf*np.maximum(0, w)

        # Total weight is product of all neutrinos in event
        unidf.columns = pd.MultiIndex.from_product([["wgt"], [colprefix], uni_columns])

        return unidf

    # All systematics together
    df = add_weights(df, make_unidf(weights.genie_systematics, "xsec"))
    df = add_weights(df, make_unidf(weights.beam_systematics, "flux"))
    df = add_weights(df, make_unidf(weights.g4_systematics, "g4"))

    # JB added 4/27/26
    # Do the cathode and gap universes by hand here
    cathode_univ = np.random.normal(size=NUNI)
    gap_univ = np.random.normal(size=NUNI)

    # Multiply together to get total
    for i in range(NUNI):
        cathode_err_wgt = np.maximum(1 + cathode_univ[i]*cathode_cross_weight_err, 0)
        gap_err_wgt = np.maximum(1 + gap_univ[i]*gap_cross_weight_err, 0)
        df[("wgt", "track_split", "univ_%i" % i, "", "", "")] = np.maximum(cathode_err_wgt*gap_err_wgt, 0)

        df[("wgt", "all", "univ_%i" % i, "", "", "")] = df[("wgt", "xsec", "univ_%i" % i, "", "", "")]*\
                                                   df[("wgt", "flux", "univ_%i" % i, "", "", "")]*\
                                                   df[("wgt", "g4", "univ_%i" % i, "", "", "")]*\
                                                   df[("wgt", "coh", "univ_%i" % i, "", "", "")]*\
                                                   df[("wgt", "track_split", "univ_%i" % i, "", "", "")]

    return Dataset(df, livetime, pot, hdrdf)

datadir = "/exp/icarus/data/users/gputnam/thesis-work/DMCP2023G/normdata/"

grl = datadir + "Run1Run2_grl.txt"
with open(grl) as f:
    goodruns = [int(l.rstrip("\n")) for l in f]

onbeam_csv = datadir + "Run1Run2_OnBeam.csv"
def onbeam_dataset(f, key, hdrkey="hdr"):
    # Data quality onbeam
    dq_cuts = [9642, 9723, 9941]
    return majority_dataset(f, key, onbeam_csv, hdrkey, dq_cuts=dq_cuts)

offbeam_csv = datadir + "Run1Run2_OffBeam.csv"
def offbeam_dataset(f, key, hdrkey="hdr"):
    # Data quality onbeam
    dq_cuts = [8518, 9781, 9837, 9851, 9892]
    return majority_dataset(f, key, offbeam_csv, hdrkey, dq_cuts=dq_cuts)

def majority_dataset(f, key, normfile, hdrkey="hdr", dq_cuts=None):
    tomrg = [
        "totpot_corr",
        "livetime_corr",
        "minbias",
        "quality_majority",
    ] + ["run", "evt"]

    normdf = pd.read_csv(normfile)
    normdf["evt"] = normdf.event
    normdf = normdf.groupby(["run", "evt"]).first().reset_index()
    hdrdf = pd.read_hdf(f, key=hdrkey)

    hdrdf = multicol_merge(hdrdf.reset_index(), normdf[tomrg], 
                on=["run", "evt"]).set_index(hdrdf.index.names)

    # Require on good run list
    hdrdf["goodrun"] = False
    for r in goodruns:
        hdrdf.loc[hdrdf.run == r, "goodrun"] = True

    # Apply data quality cut
    if dq_cuts:
        for r in dq_cuts:
            hdrdf.loc[hdrdf.run == r, "goodrun"] = False

    # cuts on events
    when_evt = hdrdf.quality_majority & ~hdrdf.minbias & hdrdf.goodrun
    # cuts on normalization
    when_norm = hdrdf.goodrun

    print(normdf.livetime_corr.sum()/1e6, hdrdf[when_norm].livetime_corr.sum()/1e6) # [s]
    print(normdf.totpot_corr.sum()/1e8, hdrdf[when_norm].totpot_corr.sum()/1e8) # POTe20

    livetime = hdrdf[when_norm].livetime_corr.sum()
    pot = hdrdf[when_norm].totpot_corr.sum()*1e12

    df = pd.read_hdf(f, key)

    return Dataset(df[broadcast(when_evt, df)], livetime, pot, hdrdf)
