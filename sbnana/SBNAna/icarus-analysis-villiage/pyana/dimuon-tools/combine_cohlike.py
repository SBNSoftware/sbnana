import argparse
import pandas as pd
import numpy as np
import sys

FNAL_COHLIKE_FILES = [
  "/icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_nom_evt.df",
  "/icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom-E_evt.df",
  "/icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom-F_evt.df",
  "/icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom-G_evt.df",
  "/icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom-H_evt.df",
]

POLARIS_COHLIKE_FILES = [
  "/icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom-polaris-prod2_evt.df",
  "/icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom-polaris-prod3_evt.df",
]

OUTPUT = "/icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom-all_nowgt_evt.df"

# Neutrinos per POT, extracted from one of the coh-like files
NU_PER_POT_BASE = 1.0448039251186124e-16

NAMES = ["mch", "evt", "hdr", "mcnu"]
with pd.HDFStore(OUTPUT) as hdf:
    for n in NAMES:
        df = None
        offset = 0
        for f in FNAL_COHLIKE_FILES: 
            thisdf = pd.read_hdf(f, n)
            if thisdf.shape[0] > 0:
                thisdf.index = thisdf.index.set_levels(thisdf.index.levels[0] + offset, level=0)
            if df is None:
                df = thisdf
            else:
                df = pd.concat([df, thisdf])
            print(f, n)
            offset += 1000000

        for f in POLARIS_COHLIKE_FILES: 
            thisdf = pd.read_hdf(f, n)
            if n == "hdr": # fix the POT for the polaris files
                thisdf.pot = 1 / NU_PER_POT_BASE
            if thisdf.shape[0] > 0:
                thisdf.index = thisdf.index.set_levels(thisdf.index.levels[0] + offset, level=0)
            df = pd.concat([df, thisdf])
            del thisdf
            print(f, n)
            offset += 1000000

        hdf.put(key=n, value=df, format="fixed")
