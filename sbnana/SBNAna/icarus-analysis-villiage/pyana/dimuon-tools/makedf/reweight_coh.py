import uproot
import numpy as np
import pandas as pd
from scipy import linalg
import util

np.random.seed(525600) # moments so dear
NUNIV = 1000

# Fit configuration file
fname = "/exp/icarus/data/users/gputnam/thesis-work/AScale0_fix.min.root"
f = uproot.open(fname)

# Load the configuration
covariance = f["covariance"].to_numpy()[0]
fit_result = f["fit_dials"].to_numpy()[0]
NVAR = fit_result.size

# Hard code "default" uncertainty
default_uncertainty = np.array([0.5, 0.9, 0.65, 0.5, 0.5, 0.5, 0.5, 1, 1])

# error on A-scaling is set separately
Aerr = 0.91393

# add in error on A to covariance
covariance = covariance + np.diag(np.concatenate(([Aerr], np.zeros((NVAR-1,)))))

# Generate the random throws for each parameter value

# The GENIE reweight convention is that the values and error are reported relative to the
# "default uncertainty" in the fit, so we need to scale by that
rands = np.random.randn(NUNIV, NVAR)
SAMPLES = 1 + (rands@linalg.cholesky(covariance) + fit_result)*default_uncertainty

# compute the A-scaling shift in each universe and the CV
A_samples = SAMPLES[:, 0]
Aexp_nom = 4/3.
Ar_A = 18
C_A = 12
Aexp_UNIV = Aexp_nom*A_samples

A_scale_UNIV = (Ar_A/C_A)**Aexp_UNIV / ((Ar_A/C_A)**Aexp_nom)

# Compute the pion scaling

# hard-code pion energy bins in x-sec reweighting
PION_E_BINS = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 6]

PION_E_SCALES = SAMPLES[:, 1:].T 
PION_E_SCALES[:4] *= A_scale_UNIV

MEAN_SCALES = np.mean(PION_E_SCALES, axis=1)

# Also add in an MA associated weight
MA_nom = 1
MA_unc = 0.3
MA_UNIV = MA_nom + MA_unc*np.random.randn(NUNIV)
def Faxial(Q2, M):
    return M**2 / (Q2 + M**2)

# And an additave uncertainty on the |t| constant
b_const_shift = 20 # 1/GeV^2
b_UNIV = b_const_shift*np.random.randn(NUNIV)

def coht(mcdf):
    nuclE = mcdf.E - mcdf.p0.genE - mcdf.p1.genE
    nuclP = mcdf.nu.P - mcdf.p0.genp - mcdf.p1.genp
    t = np.abs(nuclE**2 - util.magdf(nuclP)**2)
    return t

PION_MASS = 0.139570
def cohweight(mcdf, Q2wgt=True, bshiftwgt=True):
    wgts = []
    cccoh = (np.abs(mcdf.pdg) == 14) & (mcdf.iscc == 1) & (mcdf.genie_mode == 3)
    pionKE = mcdf.p1.genE - PION_MASS

    piwgt = pd.cut(pionKE, PION_E_BINS, labels=MEAN_SCALES).astype(float)
    wgts.append(piwgt.fillna(1).rename("cv"))
    for i in range(NUNIV):
        thiswgt = pd.cut(pionKE, PION_E_BINS, labels=PION_E_SCALES[:, i]).astype(float)
        if Q2wgt:
            thiswgt *= Faxial(mcdf.Q2, MA_UNIV[i]) / Faxial(mcdf.Q2, MA_nom) 
        if bshiftwgt:
            t = np.minimum(coht(mcdf), 0.1)
            thiswgt *= np.exp(t*b_UNIV[i])
        wgts.append(thiswgt.fillna(1).rename("univ%i" % i))

    wgts = pd.concat(wgts, axis=1)
    wgts.loc[~cccoh, :] = 1 # don't reweight non-CC-COH

    return wgts
