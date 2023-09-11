import numpy as np
import pandas as pd
from pyanalib.histogram import Histogram

def makehist(var, dataset, cut=None, bins=None, POT=None, livetime=None, areanorm=False, weight=None):
    h = Histogram(dataset, var, cut, bins, weight)
    h = h.to_pot(POT) if POT else h.to_livetime(livetime) if livetime else h.to_area() if areanorm else h
    return h

def get_systematic_variations(systs, df, hnom, histf):
    diff_his = []
    diff_los = []
    for s in systs:
        w_ps = df.wgt[(*s, "ps")]
        h_ps = histf(w_ps)

        w_ms = df.wgt[(*s, "ms")]
        h_ms = histf(w_ms)

        h_lo = np.minimum(h_ps.N, h_ms.N)
        h_hi = np.maximum(h_ps.N, h_ms.N)
        diff_hi = np.abs(h_hi - hnom.N)
        diff_lo = np.abs(hnom.N - h_lo)
        diff_his.append(diff_hi**2)
        diff_los.append(diff_lo**2)

    err_hi = np.sqrt(sum(diff_his))
    err_lo = np.sqrt(sum(diff_los))

    return err_lo, err_hi

def plotmc(fig, var, dataset, cut=None, bins=None, POT=None, areanorm=False, errorbar=False, fluxsyst=False, geniesyst=False, cvweight=True, categories=None, hist_kw={}):
    cvw = dataset.df.wgt.cv if cvweight else pd.Series(1, dataset.df.index)

    h = makehist(var, dataset, cut, bins, POT, areanorm=areanorm, weight=cvw)
    bins = h.bins

    if categories:
        hs = [makehist(var, dataset, cut & c, bins, POT, weight=cvw) for c  in categories] 
        if areanorm:
            norm = np.sum(h.N)*(h.bins[1] - h.bins[0])
            hs = [hi/norm for hi in hs] 
        fhist = fig.hist([h.centers for _ in hs], bins=h.bins, weights=[hi.N for hi in hs], **hist_kw)
    else:
        fhist = fig.hist(h.centers, bins=h.bins, weights=h.N, **hist_kw)

    if errorbar:
        fig.errorbar(h.centers, h.N, h.Nerr, linestyle="none", color="black")

    fluxsyst = fluxsyst and "fluxsyst" in dataset.df.wgt.columns.get_level_values(0)
    geniesyst = geniesyst and "geniesyst" in dataset.df.wgt.columns.get_level_values(0)
    if fluxsyst:
        systs = [("fluxsyst", f) for f in dataset.df.wgt.fluxsyst.columns.get_level_values(0)[::2]]
        histf = lambda w: makehist(var, dataset, cut, bins, POT, areanorm=areanorm, weight=cvw*w)
        ferr_lo, ferr_hi = get_systematic_variations(systs, dataset.df, h, histf)
    else: 
        ferr_lo = ferr_hi = 0
    if geniesyst:
        systs = [("geniesyst", f) for f in dataset.df.wgt.geniesyst.columns.get_level_values(0)[::2]]
        histf = lambda w: makehist(var, dataset, cut, bins, POT, areanorm=areanorm, weight=cvw*w)
        gerr_lo, gerr_hi = get_systematic_variations(systs, dataset.df, h, histf)
        gerr_lo = np.sqrt(ferr_lo**2 + gerr_lo**2)
        gerr_hi = np.sqrt(ferr_hi**2 + gerr_hi**2)
    else: 
        gerr_lo = gerr_hi = 0

    if geniesyst:
        fig.fill_between(h.bins, np.append(h.N - gerr_lo, 0), np.append(h.N + gerr_hi, 0), alpha=0.5, color="gray", step="post", hatch="x")
    if fluxsyst:
        fig.fill_between(h.bins, np.append(h.N - ferr_lo, 0), np.append(h.N + ferr_hi, 0), alpha=0.5, color="gray", step="post", hatch=".")

    return fhist 

def plotdata(fig, var, dataset, beamoff_dataset=None, cut=None, bins=None, POT=None, areanorm=False):
    h = makehist(var, dataset, cut, bins, dataset.POT)
    if beamoff_dataset:
        h = h - makehist(var, beamoff_dataset, cut, bins, livetime=dataset.livetime)

    h = h.to_pot(POT) if POT else h.to_area() if areanorm else h

    return fig.errorbar(h.centers, h.N, h.Nerr, 
        linestyle="none", color="black", marker=".")

