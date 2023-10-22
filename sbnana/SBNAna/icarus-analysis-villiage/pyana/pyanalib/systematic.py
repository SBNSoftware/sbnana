import numpy as np
import pandas as pd

class SystSpectrum(object):
    def __init__(self, hs, weights=None):
        self.hs = hs
        self.weights = weights

    def cov(self, cv):
        hcov = None
        for i in range(len(self.hs)):
            thiscov = np.outer(h.N - hs.N, h.N - hs.N)
            if hcov is None: 
                hcov = thiscov
            else:
                hcov += thiscov

        return hcov / len(self.hs)
