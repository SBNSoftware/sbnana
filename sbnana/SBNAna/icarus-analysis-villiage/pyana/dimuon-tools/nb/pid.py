import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator
import pandas as pd
import uproot
from pyanalib.variable import VAR, ARGVAR

LAr_density_gmL = 1.389875

# RECOMBINATION

# ArgoNeuT params
MODA = 0.930
MODB = 0.212
Wion = 1e3 / 4.237e7

Vps = 75058
Enom = (375./385)*Vps*1e-3 /148.25
Eshort = (368.6/378.6)*Vps*1e-3 /148.25
Efield = Enom

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

# MEAN ENEGY LOSS
mass_electron = 0.5109989461 # MeV https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
muon_mass = 105.6583745 # MeV https://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf

proton_mass = 938.272
Ival = 197.0e-6
Zval = 18.0
Aval = 39.948

Ar_molar_mass = 39.9623
Ar_ZA = 18. / Ar_molar_mass
Relec = 2.817940 * 1e-13
mole = 6.0221409*1e23
Kfactor = 4*np.pi*mole*Relec**2*mass_electron # 0.307075

def Calc_MEAN_DEDX(T, thisIval=Ival, mass=proton_mass, z=1):
    gamma = (mass+T)/mass
    beta = np.power(1.0-np.power(gamma,-2.0),0.5)
    Wmax = (2.0*mass_electron*np.power(beta,2.0)*np.power(gamma,2.0))/(1.0+2.0*gamma*(mass_electron/mass)+np.power(mass_electron/mass,2.0))

    # Medium energy 
    dens_factor = 2.0*np.log(10)*np.log10(beta*gamma)-5.2146+0.19559*np.power(3.0-np.log10(beta*gamma),3.0)
    # low energy
    dens_factor[np.log10(beta*gamma) < 0.2] = 0.
    dens_factor[beta < 1e-6] = 0.
    # high energy
    dens_factor[np.log10(beta*gamma) > 3.0] = (2.0*np.log(10)*np.log10(beta*gamma)-5.2146)[np.log10(beta*gamma) > 3.0]
    dEdx_mean = z*LAr_density_gmL*Kfactor*(Zval/Aval)*np.power(beta,-2.0)*(0.5*np.log(2.0*mass_electron*np.power(beta,2.0)*np.power(gamma,2.0)*Wmax*np.power(thisIval,-2.0))-np.power(beta,2.0)-dens_factor/2.0)

    return dEdx_mean

def Calc_RR_points(KE, dRR=0.01, mass=muon_mass, z=1, thisIval=Ival, do_recombine=False):
    thisKE = KE
    KE_points = [0.]
    RR_points = [0.]

    while thisKE > 0.0:
        deltaKE = Calc_MEAN_DEDX(np.array([thisKE]), mass=mass, z=z, thisIval=thisIval) * dRR
        RR_points.append(RR_points[-1] + dRR)
        deltaKE = deltaKE[0]
        thisKE -= deltaKE
        if do_recombine:
            deltaKE = recombination(deltaKE/dRR)*dRR
 
        KE_points.append(KE_points[-1]+deltaKE)

    # KE_points = np.array(list(reversed(KE_points[:-1])))
    KE_points = np.array(KE_points)
    KE_points = np.flip(KE_points[-1] - KE_points)
    RR_points = np.array(RR_points)

    return KE_points, RR_points

# Range to Energy
KE_points_max = 1000.
KE_points, RR_points = Calc_RR_points(KE_points_max, mass=muon_mass)
RR2KE_mu = CubicSpline(RR_points, KE_points)

KE_points, RR_points = Calc_RR_points(KE_points_max, mass=proton_mass)
RR2KE_p = CubicSpline(RR_points, KE_points)

KE_points, RR_points = Calc_RR_points(KE_points_max, mass=2*proton_mass)
RR2KE_d = CubicSpline(RR_points, KE_points)

KE_points, RR_points = Calc_RR_points(KE_points_max, mass=3*proton_mass)
RR2KE_t = CubicSpline(RR_points, KE_points)

KE_points, RR_points = Calc_RR_points(KE_points_max, mass=4*proton_mass, z=2)
RR2KE_a = CubicSpline(RR_points, KE_points)

KE_points, RR_points = Calc_RR_points(KE_points_max, mass=3*proton_mass, z=2)
RR2KE_3he = CubicSpline(RR_points, KE_points)

# Range to Charge
Q_points, RR_points = Calc_RR_points(KE_points_max, mass=proton_mass, do_recombine=True)
RR2Q_p = CubicSpline(RR_points, Q_points)

Q_points, RR_points = Calc_RR_points(KE_points_max, mass=2*proton_mass, do_recombine=True)
RR2Q_d = CubicSpline(RR_points, Q_points)

Q_points, RR_points = Calc_RR_points(KE_points_max, mass=3*proton_mass, do_recombine=True)
RR2Q_t = CubicSpline(RR_points, Q_points)

Q_points, RR_points = Calc_RR_points(KE_points_max, mass=4*proton_mass, do_recombine=True, z=2)
RR2Q_a = CubicSpline(RR_points, Q_points)

Q_points, RR_points = Calc_RR_points(KE_points_max, mass=3*proton_mass, do_recombine=True, z=2)
RR2Q_3he = CubicSpline(RR_points, Q_points)

def Calc_Q2KE_points(KE, recomb, dQ0=500, mass=muon_mass, z=1):
    thisKE = KE
    KE_points = []
    Q_points = []
    while thisKE > 0.0:
        dEdx = Calc_MEAN_DEDX(np.array([thisKE]), mass=mass, z=1)
        dQdx = recomb(dEdx)
        dx = dQ0/dQdx[0]
        deltaKE = dEdx * dx
        dQ = dQdx*dx
        Q_points.append(dQ)
        thisKE -= deltaKE[0]
        KE_points.append(thisKE)
    
    KE_points = np.flip(np.array(KE_points[:-1]), axis=0)
    Q_points = np.cumsum(np.flip(np.array(Q_points[:-1]), axis=0), axis=0)

    return KE_points, Q_points

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
def dedxdf(hitdf, gain=1., A=MODA, dqdxname="dqdx", RR_QE=5):
    df = pd.concat([dedx(gain=gain, A=A, B=hitdf.beta)(hitdf[dqdxname]), hitdf[dqdxname].rename("dqdx"), hitdf.rr, hitdf.pitch], axis=1)

    # get min/max RR
    minrr = df.groupby(level=list(range(df.index.nlevels-1))).rr.min().rename("minrr")
    maxrr = df.groupby(level=list(range(df.index.nlevels-1))).rr.max().rename("maxrr")

    df = df.join(minrr)
    df = df.join(maxrr)

    # Compute the charge-energy
    if RR_QE > 0:
        qsum = (hitdf[dqdxname]*hitdf.pitch*(hitdf.rr<RR_QE))[::-1].groupby(level=list(range(df.index.nlevels-1))).cumsum()
        Bs = hitdf.beta.groupby(level=list(range(df.index.nlevels-1))).first().values
        def this_recomb(dEdx):
            return recombination(dEdx, A=A, B=Bs)/gain
        KEs, Qs = Calc_Q2KE_points(float(RR2KE_p(RR_QE*50)), this_recomb, mass=proton_mass)

        i_func = 0
        def lookup_E(x):
            nonlocal i_func
            ind = np.searchsorted(Qs[:, i_func], x)
            E = KEs[ind]
            i_func += 1
            return E

        q_ke_p = qsum.groupby(level=list(range(df.index.nlevels-1))).transform(lookup_E)
        q_ke_p.name = "qke_p"
        df = df.join(q_ke_p)

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

@VAR
def ecal(hitdf):
    return (hitdf.dedx*hitdf.pitch).groupby(level=list(range(hitdf.index.nlevels-1))).sum()

@VAR
def erange_p(hitdf):
    rangegroup = hitdf.maxrr.groupby(level=list(range(hitdf.index.nlevels-1))).first()
    return pd.Series(RR2KE_p(rangegroup), rangegroup.index)

@ARGVAR
def erange_p_map(hitdf, RR2KE=None):
    rangegroup = hitdf.maxrr.groupby(level=list(range(hitdf.index.nlevels-1))).first()
    return pd.Series(RR2KE(rangegroup), rangegroup.index)

@ARGVAR
def ecal_p(hitdf, rrcut=2):
    Q_E = (hitdf.qke_p*(hitdf.rr < rrcut)).groupby(level=list(range(hitdf.index.nlevels-1))).max()
    R_E = (hitdf.dedx*hitdf.pitch*(hitdf.rr>=rrcut)).groupby(level=list(range(hitdf.index.nlevels-1))).sum()
    return Q_E + R_E


