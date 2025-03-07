import argparse
import pandas as pd
import numpy as np
import sys

if not len(sys.argv) > 3:
    print("Usage: python add_dfs_fixpot.py <out>.df [<in>.df]")
    sys.exit(1)

outdfs = [pd.read_hdf(f, "s") for f in sys.argv[2:]]
pots = [df.pot.iloc[0] for df in outdfs]
totpot = sum(pots)
for i in range(len(outdfs)):
    outdfs[i].weight *= pots[i] / totpot

out = pd.concat(outdfs)
out.to_hdf(sys.argv[1], key="s")
