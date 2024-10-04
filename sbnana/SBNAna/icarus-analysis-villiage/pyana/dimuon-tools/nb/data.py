import numpy as np
import pandas as pd
import math

from pyanalib.dataset import Dataset
from pyanalib.panda_helpers import *
import weights

def simple_dataset(f, key):
    df = pd.read_hdf(f, key=key)
    cv = pd.Series(1, index=df.index)
    cv.name = ("wgt", "cv")
    df = multicol_add(df, broadcast(cv, df))
    return Dataset(df, -1, -1)

def simple_dataset_df(df):
    cv = pd.Series(1, index=df.index)
    cv.name = ("wgt", "cv")
    df = multicol_add(df, broadcast(cv, df))
    return Dataset(df, -1, -1)

def mc_dataset(f, key, hdrkey="hdr", mcnukey="mcnu", syst_weights=True, mccut=None, mccut_any=True):
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

    # Central-Value
    cv = math.prod([mcdf[w] for w in weights.cv]).groupby(level=list(range(mcdf.index.nlevels-1))).prod()
    cv.name = ("wgt", "cv")
    df = multicol_add(df, broadcast(cv, df))

    if not syst_weights:
        return Dataset(df, livetime, pot, hdrdf)

    # Set the seed to be reproducible
    np.random.seed(24601)

    # Systematics
    NUNI = 100
    def make_unidf(systematics, colprefix):
        uni_columns = ["univ_%i" % i for i in range(NUNI)]
        unidf = pd.DataFrame(1, index=mcdf.index, columns=uni_columns)
        for s in systematics:
            print(s)

            # +/-1 sigma
            if mcdf[s].columns[0][0].startswith("ps"):
                ps = mcdf[s].columns[0][0]
                # ps = "ps1"
                
                shift = np.random.normal(size=NUNI)
                w = 1 + pd.DataFrame(np.outer((mcdf[s][ps]-1), shift), index=mcdf.index, columns=uni_columns)
        
            # One-sided
            elif mcdf[s].columns[0][0] == "morph":
                shift = np.random.normal(size=NUNI)
                w = 1 + pd.DataFrame(np.outer((mcdf[s].morph-1)*2, np.abs(shift)), index=mcdf.index, columns=uni_columns)                    
        
            # Universe uncertainties
            elif mcdf[s].columns[0][0] == "univ_0":
                w = mcdf[s][uni_columns]
                w.columns = [c[0] for c in w.columns]

            unidf = unidf*np.maximum(0, w)

        # Total weight is product of all neutrinos in event
        unidf = unidf.groupby(level=[0,1]).prod()
        unidf.columns = pd.MultiIndex.from_product([["wgt"], [colprefix], uni_columns])

        return unidf.groupby(level=[0,1]).prod()

    # All systematics together
    df = multicol_merge(df, make_unidf(weights.allsysts, "all"), left_index=True, right_index=True)
    df = multicol_merge(df, make_unidf(weights.genie_systematics, "xsec"), left_index=True, right_index=True)
    df = multicol_merge(df, make_unidf(weights.beam_systematics, "flux"), left_index=True, right_index=True)
    df = multicol_merge(df, make_unidf(weights.g4_systematics, "g4"), left_index=True, right_index=True)

    return Dataset(df, livetime, pot, hdrdf)

datadir = "/icarus/data/users/gputnam/DP2022P/reco/"
onbeam_csv = datadir + "Run1_Run2.csv"
def onbeam_dataset(f, key, hdrkey="hdr"):
    return majority_dataset(f, key, onbeam_csv, hdrkey)

offbeam_csv = datadir + "Run1_Run2_offbeam.csv"
def offbeam_dataset(f, key, hdrkey="hdr"):
    return majority_dataset(f, key, offbeam_csv, hdrkey)

def majority_dataset(f, key, normfile, hdrkey="hdr"):
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

    when_hdr = hdrdf.quality_majority & ~hdrdf.minbias

    livetime = hdrdf.livetime_corr[when_hdr].sum()
    pot = hdrdf.totpot_corr[when_hdr].sum()*1e12
    df = pd.read_hdf(f, key)

    return Dataset(df[broadcast(when_hdr, df)], livetime, pot, hdrdf)
