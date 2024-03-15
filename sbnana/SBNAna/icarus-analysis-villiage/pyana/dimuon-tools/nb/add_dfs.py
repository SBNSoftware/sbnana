import argparse
import pandas as pd
import numpy as np
import sys

if not len(sys.argv) > 3:
    print("Usage: python add_dfs.py <out>.df [<in>.df]")
    sys.exit(1)

out = pd.concat([pd.read_hdf(f, "s") for f in sys.argv[2:]])
out.to_hdf(sys.argv[1], key="s")
