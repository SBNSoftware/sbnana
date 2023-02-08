from pyanalib.variable import VAR, VarAccessor

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

SLCVAR = VarAccessor(["slc"])
TRKVAR = VarAccessor(["pfp", "trk"])

#def VAR(pfp):
#    return pfp.pfochar.VAR
