import numpy as np
import pandas as pd
import uproot
from pyanalib.variable import VAR, ARGVAR

# RECOMBINATION

# ArgoNeuT params
MODA = 0.930
MODB = 0.212
Wion = 1e3 / 4.237e7

Vps = 75058
Enom = (375./385)*Vps*1e-3 /148.25
Eshort = (368.6/378.6)*Vps*1e-3 /148.25
Efield = Enom
Efield = 0.4938
LAr_density_gmL = 1.389875

def recombination(dEdx, A=MODA, B=MODB, E=Efield):
    alpha = A
    beta = B / (LAr_density_gmL * E)

    dQdx = np.log(alpha + dEdx*beta) / (Wion * beta)
    return dQdx

def recombination_cor(dQdx, A=MODA, B=MODB, E=Efield):
    alpha = A
    beta = B / (LAr_density_gmL * E)
        
    dEdx = (np.exp(dQdx*Wion*beta)- alpha) / beta
        
    return dEdx

# EXPECTED dE/dx FILES
datadir = "/icarus/data/users/gputnam/"
fhist = datadir + "dEdxrestemplates.root"

profp = uproot.open(fhist)["dedx_range_pro"]
profmu = uproot.open(fhist)["dedx_range_mu"]

proton_dedx = profp.values()
proton_rr = profp.axis().edges()
proton_yerr = profp.errors(error_mode="s")
for i in range(len(proton_yerr)):
    if proton_yerr[i] < 1e-6:
        proton_yerr[i] = (proton_yerr[i-1] + proton_yerr[i+1]) / 2
    if proton_dedx[i] < 1e-6:
        proton_dedx[i] = (proton_dedx[i-1] + proton_dedx[i+1]) / 2
        
muon_dedx = profmu.values()
muon_rr = profmu.axis().edges()
muon_rr_center = profmu.axis().centers()
muon_yerr = profmu.errors(error_mode="s")

for i in range(len(muon_yerr)):
    if muon_yerr[i] < 1e-6:
        muon_yerr[i] = (muon_yerr[i-1] + muon_yerr[i+1]) / 2
    if muon_dedx[i] < 1e-6:
        muon_dedx[i] = (muon_dedx[i-1] + muon_dedx[i+1]) / 2

@ARGVAR
def scale_recombination(dedx, gain=1., A=MODA, B=MODB):
    return (recombination(dedx, A=A, B=B)/gain) / recombination(dedx)

@ARGVAR
def dedx(dqdx, gain=1., A=MODA, B=MODB):
    return recombination_cor(dqdx*gain, A, B).rename("dedx")

@ARGVAR
def dedxdf(hitdf, gain=1., A=MODA, B=MODB, dqdxname="dqdx"):
    df = pd.concat([dedx(gain=gain, A=A, B=B)(hitdf[dqdxname]), hitdf.rr], axis=1)
    minrr = df.groupby(level=list(range(df.index.nlevels-1))).rr.min().rename("minrr")
    maxrr = df.groupby(level=list(range(df.index.nlevels-1))).rr.max().rename("maxrr")

    df = df.join(minrr)
    df = df.join(maxrr)

    return df

def chi2(hitdf, exprr, expdedx, experr):
    dedx_exp = pd.cut(hitdf.rr, exprr, labels=expdedx).astype(float)
    dedx_err = pd.cut(hitdf.rr, exprr, labels=experr).astype(float)

    dedx_res = (0.04231 + 0.0001783*hitdf.dedx**2)*hitdf.dedx

    v_chi2 = (hitdf.dedx - dedx_exp)**2 / (dedx_err**2 + dedx_res**2)

    when_chi2 = (hitdf.rr < np.max(exprr)) & (hitdf.rr > hitdf.minrr) & (hitdf.rr < hitdf.maxrr) & (hitdf.dedx < 1000.)

    chi2_group = v_chi2[when_chi2].groupby(level=list(range(hitdf.index.nlevels-1)))

    return chi2_group.sum() / chi2_group.size()

@VAR
def hchi2u(hitdf):
    return chi2(hitdf, muon_rr, muon_dedx, muon_yerr)
    
@VAR
def hchi2p(hitdf):
    return chi2(hitdf, proton_rr, proton_dedx, proton_yerr)
