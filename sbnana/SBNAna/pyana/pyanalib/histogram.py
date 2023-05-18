import numpy as np
import pandas as pd

class Histogram(object):
    def __init__(self, dataset, var, cut=None, bins=None, weights=None):
        cut = cut(dataset.df) if cut else True
        v = var(dataset.df)[cut]
        weights = weights[cut] if weights is not None else None
        N, bins = np.histogram(v, bins=bins, weights=weights)
        self.N = N
        self.bins = bins
        self.centers = (bins[1:] + bins[:-1]) / 2.

        Nvar, _ = np.histogram(v, bins=bins, weights=(weights**2 if weights is not None else None))
        self.Nerr = np.sqrt(Nvar)

        self.POT = dataset.POT
        self.livetime = dataset.livetime

    def to_area(self):
        w = 1. / (np.sum(self.N)*(self.bins[1] - self.bins[0]))
        return self.scale(w)

    def to_pot(self, POT):
        w = POT / self.POT
        self.POT = POT
        return self.scale(w)

    def to_livetime(self, livetime):
        w = livetime / self.livetime
        self.livetime = livetime
        return self.scale(w)

    def scale(self, w):
        self.N = self.N*w
        self.Nerr = self.Nerr*w
        return self

    def __add__(self, othr):
        self.N += othr.N
        self.Nerr = np.sqrt(self.Nerr**2 + othr.Nerr**2)
        return self

    def __sub__(self, othr):
        self.N -= othr.N
        self.Nerr = np.sqrt(self.Nerr**2 + othr.Nerr**2)
        return self

    def __mul__(self, c):
        self.N = self.N * c
        self.Nerr = self.Nerr * c
        return self

    def __div__(self, c):
        self.N = self.N / c
        self.Nerr = self.Nerr / c
        return self


