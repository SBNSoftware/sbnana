import sys
import datetime as dt
from ntuple_glob import NTupleGlob
from makedf import make_numcdf
import numpy as np

def main(output, inputs):
    ntuples = NTupleGlob(inputs, None)
    df = ntuples.dataframe(f=make_numcdf, nproc="auto")
    df.to_hdf(output, key="df", mode="w")

if __name__ == "__main__":
    printhelp = len(sys.argv) < 3 or sys.argv[1] == "-h"
    if printhelp:
        print("Usage: python make_evtdf.py [output.df] [inputs.root,]")
    else:
        main(sys.argv[1], sys.argv[2:])
