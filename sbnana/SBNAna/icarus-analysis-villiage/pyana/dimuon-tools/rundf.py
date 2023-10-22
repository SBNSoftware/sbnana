import sys
import os
from pyanalib.ntuple_glob import NTupleGlob
from makedf.makedf import *
import pandas as pd
import warnings
from tables import NaturalNameWarning

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=NaturalNameWarning)

def main(output, inputs):
    ntuples = NTupleGlob(inputs, None)

    dfs = ntuples.dataframes(nproc="auto", fs=DFS)
    with pd.HDFStore(output) as hdf:
        for k,df in zip(NAMES, dfs):
            try:
                hdf.put(key=k, value=df, format="fixed")
            except Exception as e:
                print("Table %s failed to save, skipping. Exception: %s" % (k, str(e)))

if __name__ == "__main__":
    printhelp = len(sys.argv) < 4 or sys.argv[1] == "-h"
    if printhelp:
        print("Usage: python rundf.py [config.py] [output.df] [inputs.root,]")
    else:
        # Don't clog up the server you're running on -- let other processes take priority
        os.nice(10)


        exec(open(sys.argv[1]).read())
        if "NAMES" not in globals() or "DFS" not in globals():
            print("ERROR: config file must define <NAMES> and <DFS>")
            exit(1)
        if len(NAMES) != len(DFS): 
            print("ERROR: <NAMES> and <DFS> must have same length")
            exit(1)

        main(sys.argv[2], sys.argv[3:])
