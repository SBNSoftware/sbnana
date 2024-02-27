import numpy as np
import pandas as pd
from pyanalib.histogram import Histogram
import matplotlib

def makehist(var, dataset, cut=None, bins=None, POT=None, livetime=None, areanorm=False, weight=None):
    h = Histogram(var, dataset.POT, dataset.livetime, cut, bins, weight)
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

def plotmc(fig, var, dataset, cut=None, bins=None, POT=None, areanorm=False, errorbar=False, fluxsyst=False, geniesyst=False, cvweight=True, categories=None, wgtf=None, hist_kw={}):

    def do_wgtf(w):
        if wgtf: return wgtf(w)
        return w

    cvw = do_wgtf(dataset.df.wgt.cv if cvweight else pd.Series(1, dataset.df.index))

    v = var(dataset.df)
    vcut = cut(dataset.df) if cut is not None else cut
    h = makehist(v, dataset, vcut, bins, POT, areanorm=areanorm and categories is None, weight=cvw)
    bins = h.bins

    if categories:
        hs = [makehist(v, dataset, vcut & c(dataset.df) if vcut is not None else c(dataset.df), bins, POT, weight=cvw) for c in categories] 
        if areanorm:
            norm = np.sum(h.N)*(h.bins[1] - h.bins[0])
            hs = [hi.scale(1./norm) for hi in hs] 
        fhist = fig.hist([h.centers for _ in hs], bins=h.bins, weights=[hi.N for hi in hs], stacked=True, **hist_kw)
    else:
        fhist = fig.hist(h.centers, bins=h.bins, weights=h.N, **hist_kw)

    if errorbar:
        fig.errorbar(h.centers, h.N, h.Nerr, linestyle="none", color="black")

    fluxsyst = fluxsyst and "fluxsyst" in dataset.df.wgt.columns.get_level_values(0)
    geniesyst = geniesyst and "geniesyst" in dataset.df.wgt.columns.get_level_values(0)
    if fluxsyst:
        systs = [("fluxsyst", f) for f in dataset.df.wgt.fluxsyst.columns.get_level_values(0)[::2]]
        histf = lambda w: makehist(v, dataset, vcut, bins, POT, areanorm=areanorm, weight=cvw*do_wgtf(w))
        ferr_lo, ferr_hi = get_systematic_variations(systs, dataset.df, h, histf)
    else: 
        ferr_lo = ferr_hi = 0
    if geniesyst:
        systs = [("geniesyst", f) for f in dataset.df.wgt.geniesyst.columns.get_level_values(0)[::2]]
        histf = lambda w: makehist(v, dataset, vcut, bins, POT, areanorm=areanorm, weight=cvw*do_wgtf(w))
        gerr_lo, gerr_hi = get_systematic_variations(systs, dataset.df, h, histf)
        gerr_lo = np.sqrt(ferr_lo**2 + gerr_lo**2)
        gerr_hi = np.sqrt(ferr_hi**2 + gerr_hi**2)
    else: 
        gerr_lo = gerr_hi = 0

    gs = None
    fs = None
    if geniesyst:
        gs = fig.fill_between(h.bins, np.append(h.N - gerr_lo, 0), np.append(h.N + gerr_hi, 0), alpha=0.5, color="gray", step="post", hatch="x")
    if fluxsyst:
        fs = fig.fill_between(h.bins, np.append(h.N - ferr_lo, 0), np.append(h.N + ferr_hi, 0), alpha=0.5, color="gray", step="post", hatch=".")

    return fhist, fs, gs

def plotdata(fig, var, dataset, beamoff_dataset=None, cut=None, bins=None, POT=None, areanorm=False):
    varon = var(dataset.df)
    cuton = cut(dataset.df) if cut is not None else cut
    h = makehist(varon, dataset, cuton, bins, dataset.POT)
    if beamoff_dataset:
        varoff = var(beamoff_dataset.df)
        cutoff = cut(beamoff_dataset.df) if cut is not None else cut
        h = h - makehist(varoff, beamoff_dataset, cutoff, bins, livetime=dataset.livetime)

    h = h.to_pot(POT) if POT else h.to_area() if areanorm else h

    d = fig.errorbar(h.centers, h.N, h.Nerr, 
        linestyle="none", color="black", marker=".")

    return d

def makelegend(fig, mc, data=None, mc_labels=["MC"], split=False):
    (_, _, fh), fs, gs = mc
    if split:
        if data is not None:
            ld = fig.legend(handles=[data], labels=["Data"], frameon=False)
            fig.gca().add_artist(ld)

        l = fig.legend(handles=fh, labels=mc_labels, 
                    loc="upper left", title="MC Categories", bbox_to_anchor=(1, 1))
        if gs or fs:
            fig.gca().add_artist(l)

            handles = []
            labels = []
            if fs:
                handles.append(fs)
                labels.append("Flux")
            if gs:
                handles.append(gs)
                labels.append("Flux $+$ X-sec")
            lunc = fig.legend(handles=handles, labels=labels,
                    loc="lower left", title="Uncertainties", bbox_to_anchor=(1, 0))

    else:
        if not isinstance(fh, list):
            fh = [fh]
        handles = fh
        if data is not None:
            handles = [data] + handles
        labels = mc_labels
        if data is not None:
            labels = ["Data"] + labels
        if fs:
            handles.append(fs)
            labels.append("Flux Unc.")
        if gs:
            handles.append(gs)
            labels.append("Flux $+$ X-sec Unc.")
       
        l = fig.legend(handles=handles, labels=labels)

