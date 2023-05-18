import numpy as np
import pandas as pd

class Dataset(object):
    def __init__(self, df, livetime, POT):
        self.df = df
        self.POT = POT
        self.livetime = livetime

    def __getattr__(self, attr):
        return Dataset(getattr(self.df, attr), self.livetime, self.POT)

    def __call__(self, *args, **kw):
        return Dataset(self.df(*args, **kw), self.livetime, self.POT)

    def __getitem__(self, key):
        return Dataset(self.df[key], self.livetime, self.POT)
