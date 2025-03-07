import sys
import datetime as dt
from ntuple_glob import NTupleGlob
from makedf import make_evtdf, make_mcdf, make_numcdf, make_potdf, make_trunkhitdf, make_branchhitdf, make_crthitdf
import numpy as np
import pandas as pd

def main(output, inputs):
    ntuples = NTupleGlob(inputs, None)
    fs = [make_evtdf, make_mcdf, make_numcdf, make_potdf] #, make_trunkhitdf, make_branchhitdf] #, make_crthitdf]
    keys = ["evtdf", "mcdf", "numcdf", "hdr", "trunk_hitdf", "branch_hitdf"] #, "crtdf"]

    dfs = ntuples.dataframes(nproc="auto", fs=fs)
    with pd.HDFStore(output) as hdf:
        for k,df in zip(keys, dfs):
            hdf.put(key=k, value=df, format="fixed")

if __name__ == "__main__":
    printhelp = len(sys.argv) < 3 or sys.argv[1] == "-h"
    if printhelp:
        print("Usage: python make_alldfs.py [output.df] [inputs.root,]")
    else:
        main(sys.argv[1], sys.argv[2:])
