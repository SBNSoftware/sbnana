import numpy as np
import pandas as pd

def varhistogram(var, POT, livetime, cut=None, bins=10, weights=None):
    var = var[cut] if cut is not None else var
    weights = weights[cut] if weights is not None and cut is not None else weights
    N, bins = np.histogram(var, bins=bins, weights=weights)
    Nvar, _ = np.histogram(var, bins=bins, weights=(weights**2 if weights is not None else None))
    return Histogram(N, np.sqrt(Nvar), bins, POT, livetime)

class Histogram(object):
    def __init__(self, N, Nerr, bins, POT, livetime):
        self.N = N
        self.bins = bins
        self.centers = (bins[1:] + bins[:-1]) / 2.
        self.Nerr = Nerr
        self.POT = POT
        self.livetime = livetime

    def to_area(self):
        w = 1. / (np.sum(self.N)*(self.bins[1] - self.bins[0]))
        return self.scaled(w)

    def to_pot(self, POT):
        w = POT / self.POT
        return Histogram(self.N, self.Nerr, self.bins, POT, self.livetime).scaled(w)

    def to_livetime(self, livetime):
        w = livetime / self.livetime
        return Histogram(self.N, self.Nerr, self.bins, self.POT, livetime).scaled(w)

    def scaled(self, w):
        return Histogram(self.N*w, self.Nerr*w, self.bins, self.POT, self.livetime)

    def __add__(self, othr):
        if not isinstance(othr, Histogram):
            return Histogram(self.N + othr, self.Nerr, self.bins, self.POT, self.livetime)
        return Histogram(self.N + othr.N, np.sqrt(self.Nerr**2 + othr.Nerr**2), self.bins, self.POT, self.livetime)

    def __sub__(self, othr):
        if not isinstance(othr, Histogram):
            return Histogram(self.N - othr, self.Nerr, self.bins, self.POT, self.livetime)
        return Histogram(self.N - othr.N, np.sqrt(self.Nerr**2 + othr.Nerr**2), self.bins, self.POT, self.livetime)

    def __radd__(self, othr):
        return Histogram(self.N + othr, self.Nerr, self.bins, self.POT, self.livetime)
    def __rsub__(self, othr):
        return Histogram(self.N - othr, self.Nerr, self.bins, self.POT, self.livetime)

    def __mul__(self, othr):
        if not isinstance(othr, Histogram):
            return self.scaled(othr)
        return self.multiply(othr)

    def __truediv__(self, othr):
        if not isinstance(othr, Histogram):
            return self.scaled(1/othr)
        return self.divide(othr)

    def multiply(self, othr, propagate_error=True):
        err = (self.N*othr.N) * np.sqrt((self.Nerr/self.N)**2 + (othr.Nerr/othr.N)**2) if propagate_error else othr.N*self.Nerr
        return Histogram(self.N*othr.N, np.abs(err), self.bins, self.POT, self.livetime)

    def divide(self, othr, propagate_error=True):
        err = (self.N/othr.N) * np.sqrt((self.Nerr/self.N)**2 + (othr.Nerr/othr.N)**2) if propagate_error else (1/othr.N)*self.Nerr
        return Histogram(self.N/othr.N, np.abs(err), self.bins, self.POT, self.livetime)
