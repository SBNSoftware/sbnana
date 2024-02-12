import glob
import numpy as np
import uproot
import pandas as pd
from tqdm.auto import tqdm
import subprocess
from multiprocessing import Pool
import multiprocessing
import os
import dill
import sys
from functools import partial

CPU_COUNT = multiprocessing.cpu_count()

if CPU_COUNT == 0:
    CPU_COUNT = os.cpu_count()

class NTupleProc(object):
    def __init__(self, f=None, name="None"):
        self.name = name
        self.f = f

    def __call__(self, df):
        return self.f(df)

    def __bool__(self):
        return self.f is not None

    # make work with pickle for multiprocessing
    def __getstate__(self):
        return dill.dumps({"f": self.f, "name": self.name})

    def __setstate__(self, state):
        data = dill.loads(state)
        self.f = data["f"]
        self.name = data["name"]



def _loaddf(applyfs, inds,g):
    # fname, index, applyfs = inp
    index, fname = g
    # Convert pnfs to xroot URL's
    # if fname.startswith("/pnfs") and "scratch" not in fname:
    #     fname = fname.replace("/pnfs", "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr")
    # # fix xroot URL's
    # elif fname.startswith("xroot"):
    #     fname = fname[1:]

    madef = False

    # Flatten non-flat cafs
    # if "flat" not in fname.split("/")[-1].split("."):
    #     flatcaf = "/tmp/" + fname.split("/")[-1].split(".")[0] + ".flat.root"
    #     subprocess.run(["flatten_caf", fname, flatcaf],  stdout=subprocess.DEVNULL)
    #     fname = flatcaf
    #     madef = True
   
    try:
        f = uproot.open(fname, timeout=120)
    except (OSError, ValueError) as e:
        print("Could not open file (%s) due to exception. Skipping..." % fname) 
        print(e)
        return None
    with f:
        try:
            
            dfs = [applyf(f,inds[i]) for i,applyf in enumerate(applyfs)]
        except Exception as e:
            print("Error processing file (%s). Skipping..." % fname)
            print(e)
            return None

        # Set an index on the NTuple number to make sure we keep track of what is where
        for i in range(len(dfs)):
            if dfs[i] is not None:
                dfs[i]["__ntuple"] = index
                dfs[i].set_index("__ntuple", append=True, inplace=True)
                dfs[i] = dfs[i].reorder_levels([dfs[i].index.nlevels-1] + list(range(0, dfs[i].index.nlevels-1)))
            else:
                dfs[i] = []

    if madef:
        os.remove(fname)

    return dfs

class NTupleGlob(object):
    def __init__(self, g, branches,inds=None):
        if isinstance(g, list) and len(g) == 1:
            g = g[0]
        if isinstance(g, list):
            self.glob = g
        elif g.endswith(".list"):
            with open(g) as f:
                self.glob = [line.rstrip("\n") for line in f]
        else:
            self.glob = glob.glob(g)
        self.branches = branches
        self.inds = inds

    def dataframes(self, fs, maxfile=None, nproc=1, savemeta=False):
        if not isinstance(fs, list):
            fs = [fs]

        thisglob = self.glob 
        if maxfile:
            thisglob = thisglob[:maxfile]

        if nproc == "auto":
            nproc = min(CPU_COUNT, len(thisglob))

        ret = []
        with Pool(processes=nproc) as pool:
            for i, dfs in enumerate(tqdm(pool.imap_unordered(partial(_loaddf, fs, self.inds), enumerate(thisglob)), total=len(thisglob), unit="file", delay=5, smoothing=0.6)):
                if dfs is not None:
                    #[df.reset_index(level='__ntuple', drop=True, inplace=True) 
                    # for df in dfs if '__ntuple' in df.index.names]
                    ret.append(dfs)
        
        ret = [pd.concat([dfs[i] for dfs in ret], axis=0, ignore_index=False) for i in range(len(fs))] 
            
        if self.inds is not None:
            return ret
        
        else:
            # Fix the index So that we don't need __ntuple
            for i in range(len(ret)):
                sub_index = ret[i].index.names[2:]
                ret[i] = ret[i].reset_index()
                ret[i].entry = ret[i].groupby(["__ntuple", "entry"]).ngroup()
                ret[i].set_index(["entry"] + sub_index, inplace=True, verify_integrity=True)
                ret[i].sort_index(inplace=True)
                if not savemeta:
                    del ret[i]["__ntuple"]
                    
            return ret


