import numpy as np

from pyanalib.variable import VAR

from util import *
import var

@VAR
def stopping(e):
    return (e == 1) | (e == 2)

@VAR
def end_process(trk):
    return trk.truth.p.end_process

@VAR
def true_muon(trk):
    return np.abs(trk.truth.p.pdg) == 13

@VAR
def true_proton(trk):
    return np.abs(trk.truth.p.pdg) == 2212

@VAR
def true_trk_stop(trk):
    return TrkInFV(trk.truth.p.end)

@VAR
def trk_from_nu(trk):
    return trk.truth.p.interaction_id >= 0

@VAR
def slc_from_nu(slc):
    return slc.tmatch.idx >= 0

trk_mc_labels = ["Contained $\\nu$ $\\mu$", 
              "Exiting $\\nu$ $\\mu$",
              "Cosmic $\\mu$", 
              "Stopping $p$",
              "Other $p$",
              #"$\\pi$",
              "Other"]
    
trk_mc_categories = [
        true_muon & true_trk_stop & trk_from_nu,
        true_muon & ~true_trk_stop & trk_from_nu,
        true_muon & ~trk_from_nu,
        true_proton & true_trk_stop & stopping@end_process,
        true_proton & (~true_trk_stop | ~stopping@end_process),
        ~true_muon & ~true_proton 
]

mc_nu_cosmic_labels = ["$\\nu$", "Cosmic"]
mc_nu_cosmic_categories = [slc_from_nu, ~slc_from_nu] 

@VAR
def fiducial(slc):
    return SlcInFV(slc.vertex)

@VAR
def is_proton(trk):
    return (var.chi2p(trk) < 60) & (var.chi2u(trk) > 40)
