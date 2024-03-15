import numpy as np
import pandas as pd
import math
import uproot

from pyanalib.dataset import Dataset
from pyanalib.panda_helpers import *
import weights
import reweight_coh
import numiweight

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

def mc_dataset(f, key, hdrkey="hdr", mcnukey="mcnuwgt", syst_weights=True, mccut=None, mccut_any=True, isscalar=False):
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

    # Central-Value, if we can
    try: 
        cv = math.prod([mcdf[w] for w in weights.cv])
        cv.name = ("wgt", "cv")
        df = add_weights(df, pd.DataFrame(cv))
    except:
        if syst_weights:
            raise
        else:
            df = multicol_add(df, pd.Series(1, index=df.index, name=("wgt", "cv")))

    # load concrete correction
    nuE = mcdf.E

    # apply the concrete weight differently for nuetrinos and scalars
    # Neutrinos: use the "total" for that pdg
    # Scalars: lookup the weight for the corresponding parent pdg
    # re-weight the 
    nupdg = mcdf.parent_pdg if isscalar else mcdf.pdg

    fluxcorr_wgt = numiweight.concrete_cv(nupdg, nuE) 
    fluxcorr_wgt.name = ("wgt", "concrete", "cv", "", "", "")
    df = add_weights(df, pd.DataFrame(fluxcorr_wgt))
    df[("wgt", "cv", "", "", "", "")] *= df[("wgt", "concrete", "cv", "", "", "")]

    print("Generating Coh-weights!")
    # generate coherent pion weights
    cohweight = reweight_coh.cohweight(mcdf, douniv=syst_weights)
    cohweight.columns = pd.MultiIndex.from_tuples([("wgt", "coh", "cv", "", "", "")] + 
                                                  [("wgt", "coh", "univ_%i" % i, "", "", "") for i in range(len(cohweight.columns) - 1)])
    df = add_weights(df, cohweight)
    df[("wgt", "cv", "", "", "", "")] *= df[("wgt", "coh", "cv", "", "", "")]
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
            if s not in mcdf.columns: continue

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
        unidf.columns = pd.MultiIndex.from_product([["wgt"], [colprefix], uni_columns])

        return unidf

    # All systematics together
    df = add_weights(df, make_unidf(weights.genie_systematics, "xsec"))
    df = add_weights(df, make_unidf(weights.beam_systematics, "flux"))
    df = add_weights(df, make_unidf(weights.g4_systematics, "g4"))

    # Multiply together to get total
    for i in range(NUNI):
        df[("wgt", "all", "univ_%i" % i, "", "", "")] = df[("wgt", "xsec", "univ_%i" % i, "", "", "")]*\
                                                   df[("wgt", "flux", "univ_%i" % i, "", "", "")]*\
                                                   df[("wgt", "g4", "univ_%i" % i, "", "", "")]*\
                                                   df[("wgt", "coh", "univ_%i" % i, "", "", "")]

    return Dataset(df, livetime, pot, hdrdf)

datadir = "/icarus/data/users/gputnam/DMCP2023G/normdata/"

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
