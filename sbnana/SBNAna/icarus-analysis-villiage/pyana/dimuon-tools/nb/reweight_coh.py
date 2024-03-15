import uproot
import numpy as np
import pandas as pd
from scipy import linalg
import util

np.random.seed(525600) # moments so dear
NUNIV = 100

# Fit configuration file
# fname = "/icarus/data/users/gputnam/AScale0_fix.min.root" # pre b-scale fix
fname = "/icarus/data/users/gputnam/fitresult.min.root"
f = uproot.open(fname)

# Load the configuration
covariance = f["covariance"].to_numpy()[0]
fit_result = f["fit_dials"].to_numpy()[0]
NVAR = fit_result.size

# Hard code "default" uncertainty
default_uncertainty = np.array([0.5, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1])

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
Ar_A = 40
C_A = 12
Aexp_UNIV = Aexp_nom*A_samples

A_scale_UNIV = (Ar_A/C_A)**Aexp_UNIV / ((Ar_A/C_A)**Aexp_nom)

# BS and RS b-values
HBARC = 0.197326980 # GeV*fm
NUCLEAR_RADIUS = 1.00 # fm, in GENIE
NUCLEAR_RADIUS = NUCLEAR_RADIUS / HBARC # 1/GeV
bRS_C = (1./3)*(NUCLEAR_RADIUS*(C_A**(1./3)))**2 # Rein Segall, on Carbon
bRS = (1./3)*(NUCLEAR_RADIUS*(Ar_A**(1./3)))**2 # Rein Segall
b_BS_PION_E_BINS = [0, 0.076, 0.080, 0.100, 0.148, 0.162, 0.226, 0.486, 0.584, 0.662, 0.776, 0.870, np.inf]
b_BS_VALS_C = np.array([116, 109, 89.8, 91.0, 89.2, 80.8, 54.6, 55.2, 58.4, 60.5, 62.2, bRS_C]) # b values on carbon
b_BS_VALS = b_BS_VALS_C * (Ar_A / C_A)**(2./3) # assume A-scaling to convert to argon

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

# divide out the mean effect to get the variations
PION_E_SCALES = (PION_E_SCALES.T / MEAN_SCALES).T

def coht(mcdf):
    nuclE = mcdf.E - mcdf.p0.genE - mcdf.p1.genE
    t = np.abs(nuclE**2 - util.dmagdf(mcdf.momentum, mcdf.p0.genp+mcdf.p1.genp)**2)
    return t

PION_MASS = 0.139570
def cohweight(mcdf, Q2wgt=True, fixb=True, douniv=True):
    wgts = []
    cccoh = (np.abs(mcdf.pdg) == 14) & (mcdf.iscc == 1) & (mcdf.genie_mode == 3)
    pionKE = np.maximum(mcdf.p1.genE - PION_MASS, 0)

    piwgt = pd.cut(pionKE, PION_E_BINS, labels=MEAN_SCALES).astype(float)
    wgts.append(piwgt.fillna(1).rename("cv"))

    # GENIE by default uses the Rein-Segall b's. Correct these to Berger-Segall
    fixb_wgt = 1
    if fixb:
        t = coht(mcdf)
        bBS = pd.cut(pionKE, b_BS_PION_E_BINS, labels=b_BS_VALS).astype(float)
        fixb_wgt = np.exp(-t*(bBS - bRS))*(bBS/bRS)**2
    wgts[0] *= fixb_wgt

    if not douniv:
        wgts = pd.concat(wgts, axis=1)
        wgts.loc[~cccoh, :] = 1 # don't reweight non-CC-COH

        return wgts

    for i in range(NUNIV):
        thiswgt = np.maximum(0, pd.cut(pionKE, PION_E_BINS, labels=PION_E_SCALES[:, i]).astype(float).fillna(0))
        if Q2wgt:
            thiswgt *= Faxial(mcdf.Q2, MA_UNIV[i]) / Faxial(mcdf.Q2, MA_nom) 
        wgts.append(thiswgt.fillna(1).rename("univ%i" % i))

    wgts = pd.concat(wgts, axis=1)
    wgts.loc[~cccoh, :] = 1 # don't reweight non-CC-COH

    return wgts
