import numpy as np
import pandas as pd
from pyanalib.histogram import varhistogram
from pyanalib.variable import SystVariable
import matplotlib
import weights

class SystSpectrum(object):
    def __init__(self, hs, hscales=None):
        self.hs = hs
        self.hscales = hscales

    def cov(self, cv):
        hcov = None
        for i in range(len(self.hs)):
            thisscale = self.hscales[i] if self.hscales else 1

            thiscov = np.outer(thisscale*(cv.N - self.hs[i].N), thisscale*(cv.N - self.hs[i].N))

            if hcov is None:
                hcov = thiscov
            else:
                hcov += thiscov

        return hcov / len(self.hs)

def makehist(var, dataset, cut=None, bins=None, 
  POT=None, livetime=None, areanorm=False, 
  categories=None,
  cvweight=True, wgtf=None, 
  systematics=False, syst_weightname="all", syst_datasets=None, syst_dataset_cut=None):

    if wgtf is None:
        wgtf = lambda x: x

    cvw = dataset.df.wgt.cv if cvweight else pd.Series(1, dataset.df.index)

    if isinstance(cut, SystVariable):
        c0 = cut.cv()
    else:
        c0 = cut
    c = c0(dataset.df) if c0 else pd.Series(True, dataset.df.index)

    if isinstance(var, SystVariable):
        v0 = var.cv()
    else:
        v0 = var
    v = v0(dataset.df)
    
    def histf(w=1., do_areanorm=True, thiscut=c, thisdataset=dataset, thisvar=v, thiscv=cvw):
        h = varhistogram(thisvar, thisdataset.POT, thisdataset.livetime, thiscut, bins, wgtf(thiscv*w))
        h = h.to_pot(POT) if POT else h.to_livetime(livetime) if livetime else h.to_area() if areanorm and do_areanorm else h
        return h

    if categories:
        h = histf(do_areanorm=False)

        hs = [histf(thiscut=c&s(dataset.df), do_areanorm=False) for s in categories]
        if areanorm: # do it across all categories
            norm = np.sum(h.N)*(h.bins[1] - h.bins[0])
            hs = [hi.scaled(1./norm) for hi in hs] 

        h = hs
    else:
        h = histf()

    hsyst = []
    if systematics:
        # weights are all multisims
        if syst_weightname is not None: 
            hsyst_weights = []
            for col in dataset.df.wgt[syst_weightname].columns:
                if not col[0].startswith("univ"): continue
                hsyst_weights.append(histf(dataset.df.wgt[syst_weightname][col]))
            hsyst.append(SystSpectrum(hsyst_weights))

        # variations on the cut
        if isinstance(cut, SystVariable):
            for cprime in cut.systs():
                 cp = cprime(dataset.df)
                 if isinstance(cp, list): # handle multi-var
                     hsyst.append(SystSpectrum([histf(thiscut=cps) for cps in cp]))
                 else:
                     hsyst.append(SystSpectrum([histf(thiscut=cp)]))

        # variations on the variable
        if isinstance(var, SystVariable):
            for vprime in var.systs():
                 vp = vprime(dataset.df)
                 if isinstance(cp, list): # handle multi-var
                     hsyst.append(SystSpectrum([histf(thisvar=vps) for vps in vp]))
                 else:
                     hsyst.append(SystSpectrum([histf(thisvar=vp)]))

        # variations on the dataset
        if syst_datasets is not None:
            systc = syst_dataset_cut(dataset.df) if syst_dataset_cut else pd.Series(True, dataset.df.index)
            h_all = histf()
            h_syst_subset = histf(thiscut=c&systc)
            h_scale = h_all.divide(h_syst_subset)

            # set default weights
            for dprimes in syst_datasets:
                hds = []
                for dprime in dprimes:
                    c = c0(dprime.df) if cut else pd.Series(True, dprime.df.index)
                    v = v0(dprime.df)
                    cvw_prime = dprime.df.wgt.cv if cvweight else pd.Series(1, dprime.df.index)
                hds.append(histf(thisdataset=dprime, thisvar=v, thiscut=c, thiscv=cvw_prime))

            hsyst.append(SystSpectrum(hds, [h_scale.N for _ in hds]))

    return h, hsyst

def makecovariance(h, hsyst):
    return sum([s.cov(h) for s in hsyst])

def plotmc(fig, hs, scov=None, **hist_kw):
    if not isinstance(hs, list):
        hs = [hs]

    fhist = fig.hist([h.centers for h in hs], bins=hs[0].bins, weights=[hi.N for hi in hs], stacked=True, **hist_kw)

    if scov is not None:
        htot = sum(hs)
        fill = fig.fill_between(htot.bins, np.append(htot.N - np.sqrt(np.diag(scov)), 0), np.append(htot.N + np.sqrt(np.diag(scov)), 0), 
            alpha=0.5, color="gray", step="post", hatch="x")
    else: fill = None

    return fhist, fill

def plotdata(fig, h, **hist_kw):
    d = fig.errorbar(h.centers, h.N, h.Nerr, 
        linestyle="none", color="black", marker=".", **hist_kw)

    return d

def plotratio(fig, num, denom, scov=None):
    h = num.divide(denom, False)
    d = fig.errorbar(h.centers, h.N, np.abs(num.Nerr/denom.N), 
        linestyle="none", color="black", marker=".")

    if scov is not None:
        fill = fig.fill_between(h.bins, np.append(1 - np.sqrt(np.diag(scov))/denom.N, 0), np.append(1 + np.sqrt(np.diag(scov))/denom.N, 0), 
            alpha=0.5, color="gray", step="post", hatch="x")
    else: fill = None

    return d, fill


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

