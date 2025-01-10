from pyanalib.variable import VAR
from pyanalib.panda_helpers import *
from util import *

# Re-export stuff
from pid import dedxdf, dedx, hchi2u, hchi2p, scale_recombination

@VAR
def vabs(s):
    return np.abs(s)

@VAR
def npfp(df):
    group = list(range(df.index.nlevels-1))
    s = df.groupby(level=group).size()
    return broadcast(s, df)

@VAR
def nlongtrk(df, lencut=25):
    group = list(range(df.index.nlevels-1))
    s = (df.pfp.trk.len > lencut).groupby(level=group).sum()
    return broadcast(s, df)

@VAR
def trk(df):
    return df.pfp.trk

@VAR
def slc(df):
    return df.slc

@VAR
def pfp(df):
    return df.pfp

@VAR
def chi2p(trk):
    return trk.chi2pid.I2.chi2_proton

@VAR
def chi2u(trk):
    return trk.chi2pid.I2.chi2_muon

@VAR
def crlongtrkdiry(slc):
    return slc.nuid.crlongtrkdiry

@VAR
def trkscore(pfp):
    return pfp.trackScore

@VAR
def chgendfrac(pfp):
    return pfp.pfochar.chgendfrac

@VAR
def chgfracspread(pfp):
    return pfp.pfochar.chgfracspread

@VAR
def linfitdiff(pfp):
    return pfp.pfochar.linfitdiff

@VAR
def recop(pfp):
    ret = pfp.trk.rangeP.p_muon
    ret[~TrkInFV(pfp.trk.end)] = pfp.trk.mcsP.fwdP_muon[~TrkInFV(pfp.trk.end)]
    return ret

@VAR
def SLCVAR(df):
    return df.slc

@VAR
def TRKVAR(df):
    return df.pfp.trk

@VAR
def DF(df):
    return df
