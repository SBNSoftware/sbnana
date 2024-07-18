import numpy as np
from scipy import interpolate

# CONSTANTS

# RECOMBINATION
LAr_density_gmL_mc = 1.390

# ArgoNeuT params
MODA_mc = 0.930
MODB_mc = 0.212
Wion = 1e3 / 4.237e7
Efield_mc = 0.4938

# Measured Recombination
Ival_data = 197.0e-6
LAr_density_gmL_data = 1.392
MODA_data = 0.906
B90 = 0.203
R = 1.25
def ellipsoid_beta(phi, B90=B90, R=R):
    return B90 / np.sqrt(np.sin(phi)**2 + np.cos(phi)**2/R**2)
MODB_data = B90
Efield_data = 0.4926

# other parameters
mass_electron = 0.5109989461 # MeV https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
muon_mass = 105.6583745 # MeV https://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf

proton_mass = 938.272
Zval = 18.0
Aval = 39.948

Ival_mc = 188e-6

Ar_molar_mass = 39.9623
Ar_ZA = 18. / Ar_molar_mass
Relec = 2.817940 * 1e-13
mole = 6.0221409*1e23
Kfactor = 4*np.pi*mole*Relec**2*mass_electron # 0.307075

def recombination_cor(dQdx, A=MODA_mc, B=MODB_mc, E=Efield_mc, rho=LAr_density_gmL_mc):
    alpha = A
    beta = B / (rho * E)

    dEdx = (np.exp(dQdx*Wion*beta)- alpha) / beta

    return dEdx

def recombination(dEdx, A=MODA_mc, B=MODB_mc, E=Efield_mc, rho=LAr_density_gmL_mc):
    alpha = A
    beta = B / (rho * E)

    dQdx = np.log(alpha + dEdx*beta) / (Wion * beta)
    return dQdx

# Charge to Energy Reco
def Calc_MEAN_DEDX(T, thisIval=Ival_mc, rho=LAr_density_gmL_mc, mass=proton_mass, z=1):
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
    dEdx_mean = z*rho*Kfactor*(Zval/Aval)*np.power(beta,-2.0)*(0.5*np.log(2.0*mass_electron*np.power(beta,2.0)*np.power(gamma,2.0)*Wmax*np.power(thisIval,-2.0))-np.power(beta,2.0)-dens_factor/2.0)

    return dEdx_mean

def Calc_Q2KE_points(KE, recomb, dQ0=500, mass=proton_mass, z=1):
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

def make_interp(x, y):
    return interpolate.interp1d(x.flatten(), y.flatten(), axis=0, bounds_error=False, kind="linear", fill_value=(y[0], y[-1]))

recombination_mc = recombination
KEs, Qs = Calc_Q2KE_points(1000, recombination_mc)
Q2KE_mc = make_interp(Qs, KEs)
