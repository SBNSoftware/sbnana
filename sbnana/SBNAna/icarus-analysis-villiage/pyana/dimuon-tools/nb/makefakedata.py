import argparse
import pandas as pd
import numpy as np
import sys

output_file = sys.argv[1]
input_files = sys.argv[2:]

COHWEIGHT = np.random.random() + 0.5
SIG_SCALE = np.random.random()*50

# print(COHWEIGHT, SIG_SCALE)

outdfs = []
for i, f in enumerate(input_files):
    df = pd.read_hdf(f, "s")

    df.loc[df.genie_mode == 3, "weight"] *= COHWEIGHT
    if i == 1: # signal MC
        df.weight *= SIG_SCALE / df.weight.sum()

    draws = np.random.poisson(df.weight).astype(int)
    # print(f, df.weight.sum(), draws.sum())

    cols_2_save = [
      "mass",
      "thnumi",
      "openangle",
      "tuu"
    ]

    df = pd.DataFrame(np.repeat(df.values, draws, axis=0), columns=df.columns)[cols_2_save].astype(float)
    # df["ifile"] = i
    # print(df[df.thnumi < 5])
    outdfs.append(df)

pd.concat(outdfs).reset_index(drop=True).to_hdf(output_file, key="s", format="fixed")
