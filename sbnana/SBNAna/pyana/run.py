import sys
from pyanalib.ntuple_glob import NTupleGlob
from makedf.makedf import *
import pandas as pd
import argparse

def main(output, inputs):
    ntuples = NTupleGlob(inputs, None, INDS)
    dfs = ntuples.dataframes(nproc="auto", fs=DFS)
    with pd.HDFStore(output) as hdf:
        for k,df in zip(NAMES, dfs):
            hdf.put(key=k, value=df, format="fixed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some inputs.')
    parser.add_argument('-i', '--input', help='Input file name', required=True)
    parser.add_argument('-o', '--output', help='Output file name', required=True)
    parser.add_argument('-c', '--config', help='Config file name', required=True)
    args = parser.parse_args()

    exec(open(args.config).read())
    if "NAMES" not in globals() or "DFS" not in globals() or "INDS" not in globals():
        print("ERROR: config file must define <NAMES>, <INDS> and <DFS>")
        exit(1)
    if len(NAMES) != len(DFS):
        print("ERROR: <NAMES> and <DFS> must have same length")
        exit(1)
    main(args.output, args.input)
