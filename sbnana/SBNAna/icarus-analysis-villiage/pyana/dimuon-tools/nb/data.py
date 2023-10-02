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

def mc_dataset(f, key, hdrkey="hdr", mcnukey="mcnu"):
    hdrdf = pd.read_hdf(f, key=hdrkey)
    livetime = 0.
    pot = hdrdf.pot.sum()
    df = pd.read_hdf(f, key)

    # add in weights
    mcdf = pd.read_hdf(f, key=mcnukey)

    # Central-Value
    cv = math.prod([mcdf[w] for w in weights.cv]).groupby(level=list(range(mcdf.index.nlevels-1))).prod()
    cv.name = ("wgt", "cv")
    df = multicol_add(df, broadcast(cv, df))

    # Flux systematics
    for f in weights.fluxsyst:
        w = mcdf[f].groupby(level=list(range(mcdf.index.nlevels-1))).prod()
        w.name = tuple(["wgt", "fluxsyst"] + list(f))

        df = multicol_add(df, broadcast(w, df))

    # Genie systematics
    geniesyst = [c for c in mcdf.columns if c[0].startswith("GENIEReWeight_ICARUS_v1_multisigma")]
    for f in geniesyst:
        w = mcdf[f].groupby(level=list(range(mcdf.index.nlevels-1))).prod()
        w.name = tuple(["wgt", "geniesyst"] + list(f))
        
        df = multicol_add(df, broadcast(w, df))
    
    return Dataset(df, livetime, pot)

datadir = "/icarus/data/users/gputnam/DP2022P/reco/"
onbeam_csv = datadir + "Run1.csv"
def onbeam_dataset(f, key, hdrkey="hdr"):
    return majority_dataset(f, key, onbeam_csv, hdrkey)

offbeam_csv = datadir + "Run1_offbeam.csv"
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

    hdrdf = pd.merge(hdrdf.reset_index(), normdf[tomrg], 
                on=["run", "evt"]).set_index(hdrdf.index.names)

    when_hdr = hdrdf.quality_majority & ~hdrdf.minbias

    livetime = hdrdf.livetime_corr[when_hdr].sum()
    pot = hdrdf.totpot_corr[when_hdr].sum()*1e12
    df = pd.read_hdf(f, key)

    return Dataset(df[broadcast(when_hdr, df)], livetime, pot)
