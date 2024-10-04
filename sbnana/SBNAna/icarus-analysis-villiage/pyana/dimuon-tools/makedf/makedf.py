from pyanalib.panda_helpers import *
from .branches import *
from .util import *
from . import numisyst, g4syst, geniesyst_mcnuphase2, geniesyst_regen
import uproot
from scipy import interpolate

def make_hdrdf(f):
    hdr = loadbranches(f["recTree"], hdrbranches).rec.hdr
    return hdr

def make_mchdrdf(f):
    hdr = loadbranches(f["recTree"], mchdrbranches).rec.hdr
    return hdr

def make_potdf(f):
    pot = loadbranches(f["recTree"], potbranches).rec.hdr.numiinfo
    return pot

def make_mcnuwgtdf(f):
    return make_mcnudf(f, include_weights=True)

def make_mcnudf(f, is_regen=True, include_weights=False):
    mcdf = make_mcdf(f)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    if include_weights:
        if is_regen:
            wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst_regen.geniesyst(f, mcdf.ind), g4syst.g4syst(f, mcdf.ind)], axis=1)
        else:
            wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst_mcnuphase2.geniesyst(f, mcdf.ind)], axis=1)
        mcdf = multicol_concat(mcdf, wgtdf)
    # mcdf.index = mcdf.index.droplevel(-1)
    return mcdf

def make_mchdf(f, include_weights=False):
    mcdf = loadbranches(f["recTree"], mchbranches).rec.mc.prtl
    if include_weights:
        wgtdf = numisyst.numisyst(14, mcdf.E) # TODO: what PDG?
        mcdf = pd.concat([mcdf, wgtdf], axis=1)
    return mcdf

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

# dE/dx UNCERTAINTY FILES
func = datadir + "run1_data_dedx_uncertainty.txt"
with open(func) as f:
    philine = next(f)
    Eunc_phis = [int(p) for p in philine.rstrip("\n").split("\t") if p]
    Eunc_dedxs = []
    Euncs = []

    for line in f:
        dat = line.rstrip("\n").split("\t")
        this_dedx = float(dat[0])
        Eunc_dedxs.append(this_dedx)
        this_uncs = np.array(list(map(float, dat[1:]))) 
        this_uncs = (this_uncs/100.) * this_dedx # convert from percent to absolute
        Euncs.append(this_uncs)

dedx_phi_2dedxerr = interpolate.RegularGridInterpolator((np.array(Eunc_dedxs), np.array(Eunc_phis)), np.array(Euncs), bounds_error=False)

func = datadir + "run1_data_dedx_uncertainty_nophicorr.txt"
with open(func) as f:
    philine = next(f)
    Eunc_phis = [int(p) for p in philine.rstrip("\n").split("\t") if p]
    Eunc_dedxs = []
    Euncs = []

    for line in f:
        dat = line.rstrip("\n").split("\t")
        this_dedx = float(dat[0])
        Eunc_dedxs.append(this_dedx)
        this_uncs = np.array(list(map(float, dat[1:]))) 
        Euncs.append(this_uncs)

dedx_phi_2dedxpercenterr_nophicorr = interpolate.RegularGridInterpolator((np.array(Eunc_dedxs), np.array(Eunc_phis)), np.array(Euncs), bounds_error=False)

# RECOMBINATION
LAr_density_gmL_mc = 1.390

# ArgoNeuT params
MODA_mc = 0.930
MODB_mc = 0.212
Wion = 1e3 / 4.237e7
Efield_mc = 0.4938

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

# Data parameters
Ival_data = 197.0e-6
LAr_density_gmL_data = 1.392
MODA_data = 0.906
def ellipsoid_beta(phi, B90=0.203, R=1.25):
    return B90 / np.sqrt(np.sin(phi)**2 + np.cos(phi)**2/R**2)
MODB_data = ellipsoid_beta(1.16547) # NuMI beam direction, 66.8 degrees
Efield_data = 0.4926

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

def chi2(hitdf, exprr, expdedx, experr, dedxname="dedx"):
    dedx_exp = pd.cut(hitdf.rr, exprr, labels=expdedx).astype(float)
    dedx_err = pd.cut(hitdf.rr, exprr, labels=experr).astype(float)

    dedx_res = (0.04231 + 0.0001783*hitdf[dedxname]**2)*hitdf[dedxname]

    v_chi2 = (hitdf[dedxname] - dedx_exp)**2 / (dedx_err**2 + dedx_res**2)

    when_chi2 = (hitdf.rr < np.max(exprr)) & ~hitdf.firsthit & ~hitdf.lasthit & (hitdf[dedxname] < 1000.)

    chi2_group = v_chi2[when_chi2].groupby(level=list(range(hitdf.index.nlevels-1)))

    return chi2_group.sum() / chi2_group.size()

def chi2u(hitdf, dedxname="dedx"):
    return chi2(hitdf, muon_rr, muon_dedx, muon_yerr, dedxname)

def chi2p(hitdf, dedxname="dedx"):
    return chi2(hitdf, proton_rr, proton_dedx, proton_yerr, dedxname)

def chi2_ndof(hitdf):
    when_chi2 = (hitdf.rr < np.max(exprr)) & ~hitdf.firsthit & ~hitdf.lasthit & (hitdf.dedx < 1000.)
    chi2_group = v_chi2[when_chi2].groupby(level=list(range(hitdf.index.nlevels-1)))

    return chi2_group.size()

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

recombination_callo = lambda dEdx: recombination(dEdx)*(1 - dedx_phi_2dedxpercenterr_nophicorr((dEdx, 0))/100)
KEs, Qs = Calc_Q2KE_points(1000, recombination_callo)
Q2KE_mc_callo = make_interp(Qs, KEs)

recombination_calhi = lambda dEdx: recombination(dEdx)*(1 + dedx_phi_2dedxpercenterr_nophicorr((dEdx, 0))/100)
KEs, Qs = Calc_Q2KE_points(1000, recombination_calhi)
Q2KE_mc_calhi = make_interp(Qs, KEs)

recombination_data = lambda dEdx: recombination(dEdx, MODA_data, MODB_data, Efield_data, LAr_density_gmL_data)
KEs, Qs = Calc_Q2KE_points(1000, recombination_data) 
Q2KE_data = make_interp(Qs, KEs)

def make_trkdf(f, scoreCut=False, requiret0=False, requireCosmic=False, recalo=True):
    trkdf = loadbranches(f["recTree"], trkbranches + shwbranches)
    if scoreCut:
        trkdf = trkdf.rec.slc.reco[trkdf.rec.slc.reco.pfp.trackScore > 0.5]
    else:
        trkdf = trkdf.rec.slc.reco

    if requiret0:
        trkdf = trkdf[~np.isnan(trkdf.pfp.t0)]

    if requireCosmic:
        trkdf = trkdf[trkdf.pfp.parent == -1]

    if recalo:
        # determine MC or data
        hdrdf = make_mchdrdf(f)
        ismc = hdrdf.ismc.iloc[0]
        trkhitdf = make_trkhitdf(f)

        # compute variations if this is MC
        if ismc:
            # CALORIMETRIC VARIATIONS

            # CV_gain = 1/0.01265 # OLD 
            CV_gain = 1/0.01280 # BETTER
            CV_gain_hi = CV_gain*1.01
            CV_gain_lo = CV_gain*0.99

            CV_gain_vhi = CV_gain*1.05
            CV_gain_vlo = CV_gain*0.95

            trk_phi = np.arccos(np.abs(unitdf(trkdf.pfp.trk.truth.p.genp).x))*180/np.pi
            trk_phi = trk_phi.fillna(0) # no truth-match -- be conservative and assume most uncertain phi
            trkhitdf["phi_true"] = trk_phi
            trkdf[("pfp", "trk", "chi2pid", "phi", "", "")] = trk_phi

            trkhitdf["dedx_nom"] = recombination_cor(CV_gain*trkhitdf.dqdx*np.exp((trkhitdf.t-850)*0.4/3e3))
            trkdf[("pfp", "trk", "chi2pid", "I2", "chi2u", "")] = chi2u(trkhitdf, dedxname="dedx_nom")
            trkdf[("pfp", "trk", "chi2pid", "I2", "chi2p", "")] = chi2p(trkhitdf, dedxname="dedx_nom")

            trkhitdf["dedx_ghi"] = recombination_cor(CV_gain_hi*trkhitdf.dqdx*np.exp((trkhitdf.t-850)*0.4/3e3))
            trkdf[("pfp", "trk", "chi2pid", "I2_ghi", "chi2u", "")] = chi2u(trkhitdf, dedxname="dedx_ghi")
            trkdf[("pfp", "trk", "chi2pid", "I2_ghi", "chi2p", "")] = chi2p(trkhitdf, dedxname="dedx_ghi")

            trkhitdf["dedx_glo"] = recombination_cor(CV_gain_lo*trkhitdf.dqdx*np.exp((trkhitdf.t-850)*0.4/3e3))
            trkdf[("pfp", "trk", "chi2pid", "I2_glo", "chi2u", "")] = chi2u(trkhitdf, dedxname="dedx_glo")
            trkdf[("pfp", "trk", "chi2pid", "I2_glo", "chi2p", "")] = chi2p(trkhitdf, dedxname="dedx_glo")

            dedx_err = dedx_phi_2dedxerr(np.vstack((trkhitdf.dedx, trkhitdf.phi_true)).T)
            trkhitdf["dedx_calhi"] = trkhitdf.dedx_nom + dedx_err
            trkdf[("pfp", "trk", "chi2pid", "I2_calhi", "chi2u", "")] = chi2u(trkhitdf, dedxname="dedx_calhi")
            trkdf[("pfp", "trk", "chi2pid", "I2_calhi", "chi2p", "")] = chi2p(trkhitdf, dedxname="dedx_calhi")

            trkhitdf["dedx_callo"] = trkhitdf.dedx_nom - dedx_err
            trkdf[("pfp", "trk", "chi2pid", "I2_callo", "chi2u", "")] = chi2u(trkhitdf, dedxname="dedx_callo")
            trkdf[("pfp", "trk", "chi2pid", "I2_callo", "chi2p", "")] = chi2p(trkhitdf, dedxname="dedx_callo")

            trkhitdf["dedx_gvhi"] = recombination_cor(CV_gain_vhi*trkhitdf.dqdx*np.exp((trkhitdf.t-850)*0.4/3e3))
            trkdf[("pfp", "trk", "chi2pid", "I2_gvhi", "chi2u", "")] = chi2u(trkhitdf, dedxname="dedx_gvhi")
            trkdf[("pfp", "trk", "chi2pid", "I2_gvhi", "chi2p", "")] = chi2p(trkhitdf, dedxname="dedx_gvhi")

            trkhitdf["dedx_gvlo"] = recombination_cor(CV_gain_vlo*trkhitdf.dqdx*np.exp((trkhitdf.t-850)*0.4/3e3))
            trkdf[("pfp", "trk", "chi2pid", "I2_gvlo", "chi2u", "")] = chi2u(trkhitdf, dedxname="dedx_gvlo")
            trkdf[("pfp", "trk", "chi2pid", "I2_gvlo", "chi2p", "")] = chi2p(trkhitdf, dedxname="dedx_gvlo")

            # MCS VARIATIONS
            true_p = magdf(trkdf.pfp.trk.truth.p.genp).fillna(7.5) # be conservative, use big momentum 
            trkdf[("pfp", "trk", "mcsP_hi", "fwdP_muon", "", "")] = trkdf.pfp.trk.mcsP.fwdP_muon + true_p*1.03 # bias-hi
            trkdf[("pfp", "trk", "mcsP_lo", "fwdP_muon", "", "")] = trkdf.pfp.trk.mcsP.fwdP_muon - true_p*0.97 # bias-lo

    trkdf[("pfp", "tindex", "", "", "", "")] = trkdf.index.get_level_values(2)

    # trk_daughterdf = loadbranches(f["recTree"], pfp_daughter_branch).rec.slc.reco.pfp

    return trkdf

def make_costrkdf(f):
    trkdf = make_trkdf(f, requiret0=True, requireCosmic=True)
    slcdf = loadbranches(f["recTree"], slcbranches).rec
    return multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

def make_trkhitadcdf(f):
    return loadbranches(f["recTree"], trkhitadcbranches).rec.slc.reco.pfp.trk.calo.I2.points

def make_trkhitdf(f, include_adc=False):
    df = loadbranches(f["recTree"], trkhitbranches).rec.slc.reco.pfp.trk.calo.I2.points

    # Firsthit and Lasthit info
    ihit = df.index.get_level_values(-1)
    df["firsthit"] = ihit == 0

    lasthit = df.groupby(level=list(range(df.index.nlevels-1))).tail(1).copy()
    lasthit["lasthit"] = True
    df["lasthit"] = lasthit.lasthit
    df.lasthit = df.lasthit.fillna(False)

    # Include raw waveform information
    if include_adc:
        adcdf = make_trkhitadcdf(f)
        h_ind = adcdf.index.get_level_values(adcdf.index.nlevels-1)
        adcdf["h_t"] = broadcast(df.tdc0 - df.t, adcdf) + h_ind
        adcdf["h_w"] = broadcast(df.width, adcdf)

        # compute baselines and areas
        gl = list(range(adcdf.index.nlevels-1))
        front_baseline = adcdf.adcs[adcdf.h_t < -40].groupby(level=gl).median()
        front_baseline.name = "fped"
        back_baseline = adcdf.adcs[adcdf.h_t > 150].groupby(level=gl).median()
        back_baseline.name = "bped"

        total_area = adcdf.adcs.groupby(level=gl).sum()
        total_area.name = "tsum"
        total_count = adcdf.adcs.groupby(level=gl).size()
        total_count.name = "tnum"
        hit_area = adcdf.adcs[np.abs(adcdf.h_t) < 5*adcdf.h_w].groupby(level=gl).sum()
        hit_area.name = "hsum"
        hit_count = adcdf.adcs[np.abs(adcdf.h_t) < 5*adcdf.h_w].groupby(level=gl).size()
        hit_count.name = "hnum"

        df = multicol_add(df, front_baseline)
        df = multicol_add(df, back_baseline)
        df = multicol_add(df, total_area)
        df = multicol_add(df, total_count)
        df = multicol_add(df, hit_area)
        df = multicol_add(df, hit_count)

    return df

def make_nuclhitdf(f):
    trkhitdf = make_trkhitdf(f)
    selectdf = make_slc_trkdf(f)

    # select for nucleons
    crlongtrkdiry = selectdf.slc.nuid.crlongtrkdiry
    fiducial = SlcInFV(selectdf.slc.vertex)
    group = list(range(selectdf.index.nlevels-1))
    nlongtrk = broadcast((selectdf.pfp.trk.len > 25).groupby(level=group).sum(), selectdf)

    # Selection
    do_select = (crlongtrkdiry > -0.7) & (nlongtrk >= 1) & fiducial & (selectdf.pfp.trk.len < 10) & (selectdf.pfp.trk.chi2pid.I2.pida > 10)

    # apply onto trkhitdf
    trkhitdf = multicol_add(trkhitdf, do_select.rename("select")) 
    trkhitdf.select = trkhitdf.select.fillna(False)
    trkhitdf = trkhitdf[trkhitdf.select]
    del trkhitdf["select"]

    # Add in other helpful variables
    trkhitdf = multicol_add(trkhitdf, selectdf.pfp.trk.len.rename("len"))
    trkhitdf = multicol_add(trkhitdf, dmagdf(selectdf.pfp.trk.start, selectdf.slc.vertex).rename("dist"))
    trkhitdf = multicol_add(trkhitdf, selectdf.pfp.trk.chi2pid.I2.pida.rename("pida"))
    trkhitdf = multicol_add(trkhitdf, crlongtrkdiry.rename("crlongtrkdiry"))

    # Get timing on other planes
    for iplane in range(2):
        ts = loadbranches(f["recTree"], [trkbranch + "calo.%i.points.t" % iplane]).rec.slc.reco.pfp.trk.calo["I%i" % iplane].points.t
        tmin = ts.groupby(level=list(range(ts.index.nlevels-1))).min()
        tmax = ts.groupby(level=list(range(ts.index.nlevels-1))).max()
        trkhitdf = multicol_add(trkhitdf, tmin.rename(("P%i" % iplane, "tmin")))
        trkhitdf = multicol_add(trkhitdf, tmax.rename(("P%i" % iplane, "tmax")))

    return trkhitdf

def make_slcdf(f):
    slcdf = loadbranches(f["recTree"], slcbranches)
    slcdf = slcdf.rec

    slc_mcdf = make_mcdf(f, slc_mcbranches, slc_mcprimbranches)
    slc_mcdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "truth"] + list(c)) for c in slc_mcdf.columns])
    slcdf = multicol_merge(slcdf, slc_mcdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    return slcdf

def make_mcdf(f, branches=mcbranches, primbranches=mcprimbranches):
    # load the df
    mcdf = loadbranches(f["recTree"], branches)

    # Get the correct thing to load
    # Get the correct branches
    while mcdf.columns.nlevels > 2:
        mcdf.columns = mcdf.columns.droplevel(0)

    # Add in primary particle info
    mcprimdf = loadbranches(f["recTree"], primbranches)
    while mcprimdf.columns.nlevels > 2:
        mcprimdf.columns = mcprimdf.columns.droplevel(0)


    mcprimdf.index = mcprimdf.index.rename(mcdf.index.names[:2] + mcprimdf.index.names[2:])

    PROTON_MASS = 0.938272
    proton_KE = mcprimdf[np.abs(mcprimdf.pdg)==2212].genE.groupby(level=[0,1]).max() - PROTON_MASS
    proton_KE.name = ("max_proton_ke", "")
    mcdf = mcdf.join(proton_KE)

    mcdf.max_proton_ke = mcdf.max_proton_ke.fillna(0.)

    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==2112).groupby(level=[0,1]).sum().rename(("nn", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==2212).groupby(level=[0,1]).sum().rename(("np", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==13).groupby(level=[0,1]).sum().rename(("nmu", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==211).groupby(level=[0,1]).sum().rename(("npi", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==111).groupby(level=[0,1]).sum().rename(("npi0", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==22).groupby(level=[0,1]).sum().rename(("ng", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==321).groupby(level=[0,1]).sum().rename(("nk", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==310).groupby(level=[0,1]).sum().rename(("nk0", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==3112).groupby(level=[0,1]).sum().rename(("nsm", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==3222).groupby(level=[0,1]).sum().rename(("nsp", "")))

    # lepton info
    for b in mcprimdf.columns:
        thisb = mcprimdf[b]
        mcdf = multicol_add(mcdf, thisb.groupby(level=[0,1]).nth(0).rename(tuple(["p0"] + list(b))))
    for b in mcprimdf.columns:
        thisb = mcprimdf[b]
        mcdf = multicol_add(mcdf, thisb.groupby(level=[0,1]).nth(1).rename(tuple(["p1"] + list(b))))

    return mcdf

def make_slc_trkdf(f, trkScoreCut=False, trkDistCut=10., cutClearCosmic=True):
    # load
    trkdf = make_trkdf(f, trkScoreCut)
    slcdf = make_slcdf(f)

    # merge in tracks
    slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    # distance from vertex to track start
    slcdf = multicol_add(slcdf, dmagdf(slcdf.slc.vertex, slcdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))

    if trkDistCut > 0:
        slcdf = slcdf[slcdf.pfp.dist_to_vertex < trkDistCut]
    if cutClearCosmic:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

    return slcdf

def make_fake_evtdf(f):
    # implement the sideband selection
    def is_proton(trk):
        return (trk.chi2pid.I2.chi2_proton < 60) & (trk.chi2pid.I2.chi2_muon > 40) & TrkInFV(trk.end)
    def is_trk(pfp):
        return pfp.trackScore > 0.5
    def sideband_cut(df):
        n_proton = is_proton(df.pfp.trk).groupby(level=list(range(df.index.nlevels-1))).sum()
        trk_or_proton = is_proton(df.pfp.trk) | is_trk(df.pfp)
        n_trk_or_proton = trk_or_proton.groupby(level=list(range(df.index.nlevels-1))).sum()
        return broadcast((n_proton >= 1) & (n_trk_or_proton >= 3), df)
    # And the pre-selection
    def is_muon(df):
        trk = df.pfp.trk
        is_muon_contained = (trk.chi2pid.I2.chi2_proton > 80) & (trk.chi2pid.I2.chi2_muon < 20) & TrkInFV(trk.end)
        is_muon_exiting = ~TrkInFV(trk.end) & (trk.len > 100)
        return is_trk(df.pfp) & (is_muon_contained | is_muon_exiting)
    def evt_presel(df):
        nmuon = broadcast(is_muon(df).groupby(level=list(range(df.index.nlevels-1))).sum(), df)
        return sideband_cut(df) & (nmuon >= 2)

    evtdf = make_slc_trkdf(f)
    evtdf = evtdf[evt_presel(evtdf)].copy()

    trunk = evtdf.pfp.dist_to_vertex[is_muon(evtdf)].groupby(level=list(range(evtdf.index.nlevels-1))).idxmin()
    trunkdf = evtdf.loc[trunk].droplevel(evtdf.index.nlevels-1)

    nottrunk = evtdf.index.difference(trunk)
    branch = evtdf.pfp.dist_to_vertex.loc[nottrunk][is_muon(evtdf)].groupby(level=list(range(evtdf.index.nlevels-1))).idxmin()
    branchdf = evtdf.loc[branch].droplevel(evtdf.index.nlevels-1)

    # remove the "pfp" columns in evtdf
    for c in evtdf.columns:
        if c[0] == "pfp":
            del evtdf[c]
    # remove the "non-pfp" columns in trunk and branch
    for c in branchdf.columns:
        if c[0] != "pfp":
            del branchdf[c]
    for c in trunkdf.columns:
        if c[0] != "pfp":
            del trunkdf[c]

    # only need one entry per slice now in evtdf
    evtdf = evtdf.groupby(level=list(range(evtdf.index.nlevels-1))).first()

    # merge in the trunk and branch
    trunkdf.columns = pd.MultiIndex.from_tuples([tuple(["trunk"] + list(c[1:])) for c in trunkdf.columns])
    branchdf.columns = pd.MultiIndex.from_tuples([tuple(["branch"] + list(c[1:])) for c in branchdf.columns])

    evtdf = multicol_merge(evtdf, trunkdf, left_index=True, right_index=True, how="inner", validate="one_to_one")
    evtdf = multicol_merge(evtdf, branchdf, left_index=True, right_index=True, how="inner", validate="one_to_one")

    return evtdf

def make_stubs(f):
    stubdf = loadbranches(f["recTree"], stubbranches)
    stubdf = stubdf.rec.slc.reco.stub

    stubpdf = loadbranches(f["recTree"], stubplanebranches)
    stubpdf = stubpdf.rec.slc.reco.stub.planes

    stubdf["nplane"] = stubpdf.groupby(level=[0,1,2]).size()
    stubdf["plane"] = stubpdf.p.groupby(level=[0,1,2]).first()

    stubhitdf = loadbranches(f["recTree"], stubhitbranches)
    stubhitdf = stubhitdf.rec.slc.reco.stub.planes.hits

    stubhitdf = stubhitdf.join(stubpdf)
    stubhitdf = stubhitdf.join(stubdf.efield_vtx)
    stubhitdf = stubhitdf.join(stubdf.efield_end)

    hdrdf = make_mchdrdf(f)
    ismc = hdrdf.ismc.iloc[0]
    def dEdx2dQdx_mc(dEdx): # MC parameters
        beta = MODB_mc / (LAr_density_gmL_mc * Efield_mc)
        alpha = MODA_mc
        return np.log(alpha + dEdx*beta) / (Wion*beta)
    def dEdx2dQdx_data(dEdx): # data parameters
        beta = MODB_data / (LAr_density_gmL_data * Efield_data)
        alpha = MODA_data
        return np.log(alpha + dEdx*beta) / (Wion*beta)

    dEdx2dQdx = dEdx2dQdx_mc if ismc else dEdx2dQdx_data
    MIP_dqdx = dEdx2dQdx(1.7) 

    stub_end_charge = stubhitdf.charge[stubhitdf.wire == stubhitdf.hit_w].groupby(level=[0,1,2,3]).first().groupby(level=[0,1,2]).first()
    stub_end_charge.name = ("endp_charge", "", "")

    stub_pitch = stubpdf.pitch.groupby(level=[0,1,2]).first()
    stub_pitch.name = ("pitch", "", "")

    stubdir_is_pos = (stubhitdf.hit_w - stubhitdf.vtx_w) > 0.
    when_sum = ((stubhitdf.wire > stubhitdf.vtx_w) == stubdir_is_pos) & (((stubhitdf.wire < stubhitdf.hit_w) == stubdir_is_pos) | (stubhitdf.wire == stubhitdf.hit_w)) 
    stubcharge = (stubhitdf.charge[when_sum]).groupby(level=[0,1,2,3]).sum().groupby(level=[0,1,2]).first()
    stubcharge.name = ("charge", "", "")

    stubinccharge = (stubhitdf.charge).groupby(level=[0,1,2,3]).sum().groupby(level=[0,1,2]).first()
    stubinccharge.name = ("inc_charge", "", "")

    hit_before_start = ((stubhitdf.wire < stubhitdf.vtx_w) == stubdir_is_pos)
    stub_inc_sub_charge = (stubhitdf.charge - MIP_dqdx*stubhitdf.ontrack*(~hit_before_start)*stubhitdf.trkpitch).groupby(level=[0,1,2,3]).sum().groupby(level=[0,1,2]).first()
    stub_inc_sub_charge.name = ("inc_sub_charge", "", "")

    stubdf = stubdf.join(stubcharge)
    stubdf = stubdf.join(stubinccharge)
    stubdf = stubdf.join(stub_inc_sub_charge)
    stubdf = stubdf.join(stub_end_charge)
    stubdf = stubdf.join(stub_pitch)
    stubdf["length"] = magdf(stubdf.vtx - stubdf.end)
    stubdf["Q"] = stubdf.inc_sub_charge

    # convert charge to energy
    if ismc:
        stubdf["ke"] = Q2KE_mc(stubdf.Q)
        # also do calorimetric variations
        stubdf["ke_callo"] = Q2KE_mc_callo(stubdf.Q)
        stubdf["ke_calhi"] = Q2KE_mc_calhi(stubdf.Q)
    else:
        stubdf["ke"] = Q2KE_data(stubdf.Q)
        stubdf["ke_callo"] = np.nan
        stubdf["ke_calhi"] = np.nan

    stubdf.ke = stubdf.ke.fillna(0)
    stubdf.Q = stubdf.Q.fillna(0)

    stubdf["dedx"] = stubdf.ke / stubdf.length
    stubdf["dedx_callo"] = stubdf.ke_callo / stubdf.length
    stubdf["dedx_calhi"] = stubdf.ke_calhi / stubdf.length

    # only take collection plane
    stubdf = stubdf[stubdf.plane == 2]

    stub_length_bins = [0, 0.5, 1, 2, 3]
    stub_length_name = ["l0_5cm", "l1cm", "l2cm", "l3cm"]
    tosave = ["dedx", "dedx_callo", "dedx_calhi", "Q", "length", "charge", "inc_charge"] 

    df_tosave = []
    for blo, bhi, name in zip(stub_length_bins[:-1], stub_length_bins[1:], stub_length_name):
        stub_tosave = stubdf.dedx[(stubdf.length > blo) & (stubdf.length < bhi)].groupby(level=[0,1]).idxmax()
        for col in tosave:
            s = stubdf.loc[stub_tosave, col]
            s.name = ("stub", name, col, "", "", "")
            s.index = s.index.droplevel(-1)
            df_tosave.append(s)

    return pd.concat(df_tosave, axis=1)

def make_evt_trkhitdf(f):
    # implement the sideband selection
    def is_muon_contained(trk):
        return (trk.chi2pid.I2.chi2_proton > 80) & (trk.chi2pid.I2.chi2_muon < 30) & TrkInFV(trk.end) & (trk.len > 10)
    def is_muon_exiting(trk):
        return (trk.len > 100) & ~TrkInFV(trk.end)
    def is_muon(pfp):
        return (is_muon_contained(pfp.trk) | is_muon_exiting(pfp.trk)) & (pfp.dist_to_vertex < 10)

    trkhitdf = make_trkhitdf(f)

    # Save hits on the trunk and branch, only
    slcdf = make_slc_trkdf(f)
    trklevel = list(range(slcdf.index.nlevels-1))
    trunk = slcdf.pfp.dist_to_vertex[is_muon(slcdf.pfp)].groupby(level=trklevel).idxmin()
    nottrunk = slcdf.index.difference(trunk)
    branch = slcdf.pfp.dist_to_vertex.loc[nottrunk][is_muon(slcdf.pfp)].groupby(level=trklevel).idxmin()
    slcdf["trunk_or_branch"] = False
    slcdf.loc[trunk, "trunk_or_branch"] = True
    slcdf.loc[branch, "trunk_or_branch"] = True

    # Apply cut
    trkhitdf = multicol_add(trkhitdf, slcdf.trunk_or_branch.rename("select")) 
    trkhitdf.select = trkhitdf.select.fillna(False)
    trkhitdf = trkhitdf[trkhitdf.select]
    del trkhitdf["select"]

    return trkhitdf

def make_evtdf(f, load_hits=False):
    # implement the sideband selection
    def is_muon_contained(trk):
        return (trk.chi2pid.I2.chi2_proton > 80) & (trk.chi2pid.I2.chi2_muon < 30) & TrkInFV(trk.end) & (trk.len > 10)
    def is_muon_exiting(trk):
        return (trk.len > 100) & ~TrkInFV(trk.end)
    def is_muon(pfp):
        return (is_muon_contained(pfp.trk) | is_muon_exiting(pfp.trk)) & (pfp.dist_to_vertex < 10)
    def evt_presel(df):
        trklevel = list(range(df.index.nlevels-1))
        nmuon = broadcast(is_muon(df.pfp).groupby(level=trklevel).sum(), df)
        return SlcInFV(df.slc.vertex) & (nmuon >= 2)

    trk_slcdf = make_slc_trkdf(f, trkDistCut=-1)
    trk_slcdf = trk_slcdf[evt_presel(trk_slcdf)]
    slcdf = trk_slcdf.copy()

    trklevel = list(range(slcdf.index.nlevels-1))

    trunk = slcdf.pfp.dist_to_vertex[is_muon(slcdf.pfp)].groupby(level=trklevel).idxmin()
    trunkdf = slcdf.loc[trunk].droplevel(slcdf.index.nlevels-1)

    nottrunk = slcdf.index.difference(trunk)
    branch = slcdf.pfp.dist_to_vertex.loc[nottrunk][is_muon(slcdf.pfp)].groupby(level=trklevel).idxmin()
    branchdf = slcdf.loc[branch].droplevel(slcdf.index.nlevels-1)

    # remove the "pfp" columns in slcdf
    for c in slcdf.columns:
        if c[0] == "pfp":
            del slcdf[c]
    # remove the "non-pfp" columns in trunk and branch
    for c in branchdf.columns:
        if c[0] != "pfp":
            del branchdf[c]
    for c in trunkdf.columns:
        if c[0] != "pfp":
            del trunkdf[c]

    # only need one entry per slice now in evtdf
    evtdf = slcdf.groupby(level=trklevel).first()

    # merge in the trunk and branch
    trunkdf.columns = pd.MultiIndex.from_tuples([tuple(["trunk"] + list(c[1:])) for c in trunkdf.columns])
    branchdf.columns = pd.MultiIndex.from_tuples([tuple(["branch"] + list(c[1:])) for c in branchdf.columns])

    evtdf = multicol_merge(evtdf, trunkdf, left_index=True, right_index=True, how="inner", validate="one_to_one")
    evtdf = multicol_merge(evtdf, branchdf, left_index=True, right_index=True, how="inner", validate="one_to_one")

    # Hit information in the trunk
    if load_hits:
        trkhitdf = loadbranches(f["recTree"], trkhitbranches)
        trkhitdf = trkhitdf.rec.slc.reco.pfp.trk.calo.I2.points
        trkhitlevel = list(range(trkhitdf.index.nlevels-1))

        trunk_min_rr = evtdf.trunk.trk.len - dmagdf(evtdf.trunk.trk.start, evtdf.branch.trk.start)
        trunk_min_rr.name = ("minrr", "")
        trkhitdf = trkhitdf.join(trunk_min_rr)

        trunk_start_dqdx = trkhitdf[(trkhitdf.rr > trkhitdf.minrr) & (trkhitdf.rr > 5)].groupby(level=trkhitlevel).dqdx.median()
        trunk_start_dqdx.name = ("trunk", "med_start_dqdx", "", "", "", "")

        trunk_start_dedx = trkhitdf[(trkhitdf.rr > trkhitdf.minrr) & (trkhitdf.rr > 5)].groupby(level=trkhitlevel).dedx.median()
        trunk_start_dedx.name = ("trunk", "med_start_dedx", "", "", "", "")

        evtdf = evtdf.join(trunk_start_dqdx, how="left",
                           on=["entry", "rec.slc..index", ("trunk", "tindex")])
        evtdf = evtdf.join(trunk_start_dedx, how="left",
                           on=["entry", "rec.slc..index", ("trunk", "tindex")])

    # Get information on other objects in the slice
    not_trunk_branch = nottrunk.difference(branch)
    evtdf["third_trk_dist"] = trk_slcdf.pfp.dist_to_vertex.loc[not_trunk_branch].groupby(level=trklevel).min()
    # get the third track
    thirdtrk = trk_slcdf.pfp.dist_to_vertex.loc[not_trunk_branch].groupby(level=trklevel).idxmin().dropna()
    thirdtrkdf = trk_slcdf.loc[thirdtrk].droplevel(len(trklevel)).pfp.trk
    # evtdf["third_trk_chi2_proton"] = thirdtrkdf.chi2pid.I2.chi2_proton[thirdtrkdf.chi2pid.I2.pid_ndof > 0]
    # evtdf["third_trk_chi2_muon"] = thirdtrkdf.chi2pid.I2.chi2_muon[thirdtrkdf.chi2pid.I2.pid_ndof > 0]
    evtdf["third_trk_contained"] = TrkInFV(thirdtrkdf.end)

    # Get the other pfp df
    othrdf = trk_slcdf.loc[not_trunk_branch]
    shwdf = othrdf[othrdf.pfp.trackScore < 0.5]
    trkdf = othrdf[othrdf.pfp.trackScore >= 0.5]

    # Other-objects -- close to vertex
    chi2_ok = othrdf.pfp.trk.chi2pid.I2.pid_ndof > 0
    # min/max chi2s
    evtdf["min_othr_chi2_proton"] = othrdf[(othrdf.pfp.dist_to_vertex < 10) & chi2_ok].pfp.trk.chi2pid.I2.chi2_proton.groupby(level=trklevel).min()
    evtdf["max_othr_chi2_muon"] = othrdf[(othrdf.pfp.dist_to_vertex < 10) & chi2_ok].pfp.trk.chi2pid.I2.chi2_muon.groupby(level=trklevel).max()

    # lengths
    evtdf["max_shw_len"] = shwdf[shwdf.pfp.dist_to_vertex < 10].pfp.shw.len.groupby(level=trklevel).max()
    evtdf["max_othr_trk_len"] = trkdf[trkdf.pfp.dist_to_vertex < 10].pfp.trk.len.groupby(level=trklevel).max()

    # Other objects -- primary
    # min/max chi2s
    evtdf["min_othr_chi2_proton_primary"] = othrdf[(othrdf.pfp.parent_is_primary == 1) & chi2_ok].pfp.trk.chi2pid.I2.chi2_proton.groupby(level=trklevel).min()
    evtdf["max_othr_chi2_muon_primary"] = othrdf[(othrdf.pfp.parent_is_primary == 1) & chi2_ok].pfp.trk.chi2pid.I2.chi2_muon.groupby(level=trklevel).max()

    # lengths
    evtdf["max_shw_len_primary"] = shwdf[shwdf.pfp.parent_is_primary == 1].pfp.shw.len.groupby(level=trklevel).max()
    evtdf["max_othr_trk_len_primary"] = trkdf[trkdf.pfp.parent_is_primary == 1].pfp.trk.len.groupby(level=trklevel).max()

    # Other objects -- all
    # min/max chi2s
    evtdf["min_othr_chi2_proton_all"] = othrdf[chi2_ok].pfp.trk.chi2pid.I2.chi2_proton.groupby(level=trklevel).min()
    evtdf["max_othr_chi2_muon_all"] = othrdf[chi2_ok].pfp.trk.chi2pid.I2.chi2_muon.groupby(level=trklevel).max()

    # lengths
    evtdf["max_shw_len_all"] = shwdf.pfp.shw.len.groupby(level=trklevel).max()
    evtdf["max_othr_trk_len_all"] = trkdf.pfp.trk.len.groupby(level=trklevel).max()

    # stubs!
    evtdf = evtdf.join(make_stubs(f))

    return evtdf

    
    
