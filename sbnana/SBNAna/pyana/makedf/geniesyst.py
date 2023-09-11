import uproot
import numpy as np
import pandas as pd
from .ttree import *

def geniesyst(f, nuind):
    if "globalTree" not in f:
        return pd.DataFrame(index=nuind.index)
    treename = get_tree_name(f,'globalTree')
    nuidx = pd.MultiIndex.from_arrays([nuind.index.get_level_values(0), nuind])

    globalTree = f[treename]
    geniewgt_names = [n for n in f[treename]['global/wgts/wgts.name'].arrays(library="np")['wgts.name'][0]]
    geniewgt_types = f[treename]['global/wgts/wgts.type'].arrays(library="np")['wgts.type'][0]
    nuniv = f[treename]['global/wgts/wgts.nuniv'].arrays(library="pd").sum().values[0]

    wgts = f["recTree"]['rec.mc.nu.wgt.univ'].arrays(library="pd")
    wgts["inu"] = wgts.index.get_level_values(1) // nuniv
    wgts["iwgt"] = wgts.index.get_level_values(1) % nuniv
    wgts = wgts.reset_index().set_index(["entry", "inu", "iwgt"]).drop(columns="subentry")
    wgts.columns = ["wgt"]

    systs = []
    for i, (n, t) in enumerate(zip(geniewgt_names, geniewgt_types)):
        if t != 3: continue # only take unisim
        
        s_ps = wgts.groupby(level=[0,1]).nth(i*6).wgt
        s_ps.name = (n, "ps")

        s_ms = wgts.groupby(level=[0,1]).nth(i*6+1).wgt
        s_ms.name = (n, "ms")

        systs.append(s_ps)
        systs.append(s_ms)

    systs = pd.DataFrame(systs).T

    s_idx = systs.index.get_indexer(nuidx)
    systs_match = systs.iloc[s_idx]
    systs_match.loc[s_idx < 0, :] = 1.
    systs_match.index = nuind.index

    return systs_match

