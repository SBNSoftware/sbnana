import sys
import datetime as dt
from ntuple_glob import NTupleGlob
from makedf import make_numcdf
import numpy as np
from panda_helpers import *

def make_potdf(f):
    hdr = loadbranches(f["recTree"], ["rec.hdr.ngenevt", "rec.hdr.pot", "rec.hdr.first_in_subrun"]).rec.hdr
    return hdr

def main(output, inputs):
    ntuples = NTupleGlob(inputs, None)
    df = ntuples.dataframe(f=make_potdf, nproc="auto")
    df.to_hdf(output, key="df", mode="w")

if __name__ == "__main__":
    printhelp = len(sys.argv) < 3 or sys.argv[1] == "-h"
    if printhelp:
        print("Usage: python make_evtdf.py [output.df] [inputs.root,]")
    else:
        main(sys.argv[1], sys.argv[2:])
