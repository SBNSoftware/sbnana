import numpy as np
import pandas as pd

def reweight(time):
    goodtimes = (2, 3.5)
    badtimes = (9.5, np.inf)

    is_goodtime = (time > goodtimes[0]) & (time < goodtimes[1])
    is_badtime = (time > badtimes[0]) & (time < badtimes[1])

    timewgt_v = (is_goodtime.sum() + is_badtime.sum()) / is_goodtime.sum()
    timewgt = pd.Series(1, time.index)
    timewgt[is_goodtime] = timewgt_v
    timewgt[is_badtime] = 0
    return timewgt
