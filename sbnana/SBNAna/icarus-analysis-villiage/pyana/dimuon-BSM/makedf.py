import numpy as np
import pandas as pd
from panda_helpers import *
from scipy.interpolate import CubicSpline
import warnings
from tables import NaturalNameWarning
import numisyst
import geniesyst

SIDEBAND = False

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
warnings.simplefilter(action='ignore', category=NaturalNameWarning)

def InFV(df, inzback, inx=10, iny=10, inzfront=10):
    xmin_C0 = -358.49
    xmax_C0 = -61.94
    
    xmax_C1 = -xmin_C0
    xmin_C1 = -xmax_C0
    
    ymin = -181.85999999999999
    ymax = 134.96
    
    zmin = -894.950652270838
    zmax = 894.950652270838
    
    xmin_C0 = xmin_C0 + inx
    xmax_C0 = xmax_C0 - inx
    xmin_C1 = xmin_C1 + inx
    xmax_C1 = xmax_C1 - inx
        
    ymin = ymin + iny
    ymax = ymax - iny
    
    zmin = zmin + inzfront
    zmax = zmax - inzback
    
    return (((df.x < xmax_C0) & (df.x > xmin_C0)) | ((df.x < xmax_C1) & (df.x > xmin_C1))) &\
        (df.y < ymax) & (df.y > ymin) & (df.z < zmax) & (df.z > zmin)

def TrkInFV(df):
    return InFV(df, 15.)

def SlcInFV(df):
    return InFV(df, 100.)

def recop(trk):
    p = trk.rangeP.p_muon
    p[~TrkInFV(trk.end)] = trk.mcsP.fwdP_muon
    return p

def mag(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def magdf(df):
    return mag(df.x, df.y, df.z)
    
def dist(df1, df2):
    # return magdf(df1 - df2) 
    return mag(df1.x - df2.x, df1.y - df2.y, df1.z - df2.z)

def dca(p0, n, p1):
    return mag((p0.y - p1.y)*n.z - (p0.z - p1.z)*n.y, 
               (p0.z - p1.z)*n.x - (p0.x - p1.x)*n.z,
               (p0.x - p1.x)*n.y - (p0.y - p1.y)*n.x)

def unit(df):
    return df.divide(magdf(df), axis=0)

def dot(df1, df2):
    return df1.x*df2.x + df1.y*df2.y + df1.z*df2.z

# Compute range momentum from length
muon_range = np.array([
      9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1,
      1.063E2,  1.725E2, 2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3,
      2.297E3,  4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
      4.326E4,  5.768E4, 7.734E4, 1.060E5, 1.307E5
]) / 1.396
muon_KE = np.array([
    10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
    400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
    20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000
])
MUON_MASS = 105.6583745 # MeV

muon_range_to_KE = CubicSpline(muon_range, muon_KE)
muon_range_to_P = lambda R: np.sqrt(muon_range_to_KE(R)**2 + 2*MUON_MASS*muon_range_to_KE(R))/1e3
muon_KE_to_range = CubicSpline(muon_KE, muon_range)
muon_P_to_range = lambda P: muon_KE_to_range(np.sqrt((P*1e3)**2 + MUON_MASS**2) - MUON_MASS)

mcbranches = [
    "rec.mc.prtl.E",
    "rec.mc.prtl.position.x",
    "rec.mc.prtl.position.y",
    "rec.mc.prtl.position.z",
]

mcnubranches = [
    "rec.mc.nu.E",
    "rec.mc.nu.position.x",
    "rec.mc.nu.position.y",
    "rec.mc.nu.position.z",
    "rec.mc.nu.pdg",
    "rec.mc.nu.iscc",
    "rec.mc.nu.genie_mode"
]

mcprimbranches = [
    "rec.mc.nu.prim.genE",
    "rec.mc.nu.prim.pdg",
    "rec.mc.nu.prim.genp.x",
    "rec.mc.nu.prim.genp.y",
    "rec.mc.nu.prim.genp.z",
]


trueparticlenames = [
    "start_process",
    "end_process",
    "pdg",
    "startE",
    "start.x", "start.y", "start.z",
    "end.x", "end.y", "end.z",
    "genp.x", "genp.y", "genp.z",
    "length",
    "G4ID",
    "cont_tpc",
    "genE",
    "interaction_id"
]

trkbranch = "rec.slc.reco.pfp.trk."
shwbranch = "rec.slc.reco.pfp.shw."

pfpbranches = [
    "rec.slc.reco.pfp.parent_is_primary",
    "rec.slc.reco.pfp.slcID",
    "rec.slc.reco.pfp.trackScore",
    "rec.slc.reco.pfp.parent",
    "rec.slc.reco.pfp.id",
]

trkbranches = [
    trkbranch + "producer",
    trkbranch + "start.x", trkbranch + "start.y", trkbranch + "start.z",
    trkbranch + "end.x", trkbranch + "end.y", trkbranch + "end.z",
    trkbranch + "dir.x", trkbranch + "dir.y", trkbranch + "dir.z",
    trkbranch + "len",
    trkbranch + "rangeP.p_muon",
    trkbranch + "mcsP.fwdP_muon",
    trkbranch + "rangeP.p_pion",
    trkbranch + "mcsP.fwdP_pion",
    trkbranch + "bestplane",
    trkbranch + "crthit.distance",
    trkbranch + "crthit.hit.time",
    trkbranch + "crthit.hit.pe",
] + pfpbranches

trk_mc_branches = [
    trkbranch + "chi2pid.2.chi2_muon",
    trkbranch + "chi2pid.2.chi2_proton",
    trkbranch + "truth.pur",
    trkbranch + "truth.eff",
]

trkkinkbranches = [
    trkbranch + "kink.p.%i.cand.radius",
    trkbranch + "kink.p.%i.cand.cosopen",
    trkbranch + "kink.p.%i.cand.x",
    trkbranch + "kink.p.%i.cand.y",
    trkbranch + "kink.p.%i.cand.z",
    trkbranch + "kink.p.%i.cand.w",
]

for n in trueparticlenames: trkbranches.append(trkbranch + "truth.p." + n)

trkmatchbranches = [
    trkbranch + "truth.matches.energy"
]

trkhitbranches = [
    trkbranch + "calo.2.points.dedx",
    trkbranch + "calo.2.points.dqdx",
    trkbranch + "calo.2.points.pitch",
    trkbranch + "calo.2.points.integral",
    trkbranch + "calo.2.points.rr",
    trkbranch + "calo.2.points.wire",
    trkbranch + "calo.2.points.tpc",
    trkbranch + "calo.2.points.width",
    trkbranch + "calo.2.points.sumadc",
    trkbranch + "calo.2.points.mult",
    trkbranch + "calo.2.points.t",
    trkbranch + "calo.2.points.x",
    trkbranch + "calo.2.points.y",
    trkbranch + "calo.2.points.z",
    trkbranch + "calo.2.points.truth.e",
    trkbranch + "calo.2.points.truth.nelec",
    trkbranch + "calo.2.points.truth.pitch",
    trkbranch + "calo.2.points.truth.rr",
]

shwbranches = [
  shwbranch + "len",
  shwbranch + "start.x",
  shwbranch + "start.y",
  shwbranch + "start.z",
] + pfpbranches

slcbranches = [
    "rec.slc.charge",
    "rec.slc.is_clear_cosmic",
    "rec.slc.vertex.x", "rec.slc.vertex.y", "rec.slc.vertex.z",
    "rec.slc.self",
    "rec.slc.tmatch.eff",
    "rec.slc.tmatch.pur",
    "rec.slc.tmatch.index",
    "rec.slc.producer",
    "rec.slc.nu_score",
    "rec.slc.fmatch.score",
    "rec.slc.fmatch.time",
    "rec.slc.nuid.crlongtrkdiry"
]

slc_mcbranches = ["rec.slc.truth." + ".".join(s.split(".")[3:]) for s in mcnubranches]
slc_mcprimbranches = ["rec.slc.truth." + ".".join(s.split(".")[3:]) for s in mcprimbranches]

opflashbranches = [
    "rec.opflashes.totalpe",
    "rec.opflashes.time",
    "rec.opflashes.cryo",
    "rec.opflashes.center.y",
    "rec.opflashes.center.z",
]

hdrbranches = [
    "rec.hdr.pot",
    "rec.hdr.first_in_subrun",
    "rec.hdr.run",
    "rec.hdr.subrun",
    "rec.hdr.ngenevt",
    "rec.hdr.evt",
    "rec.hdr.proc",
    "rec.hdr.cluster",
    "rec.hdr.triggerinfo.trigger_within_gate"
]

stubbranches = [
    "rec.slc.reco.stub.vtx.x",
    "rec.slc.reco.stub.vtx.y",
    "rec.slc.reco.stub.vtx.z",
    "rec.slc.reco.stub.end.x",
    "rec.slc.reco.stub.end.y",
    "rec.slc.reco.stub.end.z",

    "rec.slc.reco.stub.efield_vtx",
    "rec.slc.reco.stub.efield_end",

    
    "rec.slc.reco.stub.truth.p.pdg",
    "rec.slc.reco.stub.truth.p.genE",
    "rec.slc.reco.stub.truth.p.interaction_id",
]

stubplanebranches = [
    "rec.slc.reco.stub.planes.p",
    "rec.slc.reco.stub.planes.hit_w",
    "rec.slc.reco.stub.planes.vtx_w",    
    "rec.slc.reco.stub.planes.pitch",    
    "rec.slc.reco.stub.planes.trkpitch",    
]

stubhitbranches = [
    "rec.slc.reco.stub.planes.hits.charge",
    "rec.slc.reco.stub.planes.hits.ontrack",
    "rec.slc.reco.stub.planes.hits.wire",
]

crthitbranches = [
    "rec.crt_hits.t0",
    "rec.crt_hits.pe",
    "rec.crt_hits.plane",
]

def make_crthitdf(f):
    crtdf = loadbranches(f["recTree"], crthitbranches).rec.crt_hits
    return crtdf[(crtdf.t0 > -20) & (crtdf.t0 < 20)]

def make_crtvetodf(f):
    crtdf = loadbranches(f["recTree"], crthitbranches).rec.crt_hits

    def tcut(t):
        return (t > 0) & (t < 10)
    def side(p):
        return (p >= 40) & (p <= 47)
    def top(p):
        return (p >= 30) & (p <= 34)

    topcrt_pe = crtdf.pe[tcut(crtdf.t0) & top(crtdf.plane)].groupby(level=0).max()
    topcrt_pe.name = "topcrt_pe"
    sidecrt_pe = crtdf.pe[tcut(crtdf.t0) & side(crtdf.plane)].groupby(level=0).max()
    sidecrt_pe.name = "sidecrt_pe"

    return pd.concat([topcrt_pe, sidecrt_pe], axis=1)

def make_potdf(f):
    hdr = loadbranches(f["recTree"], hdrbranches).rec.hdr
    return hdr

def make_flashdf(f):
    return add_header(f, loadbranches(f["recTree"], opflashbranches).rec.opflashes)

def make_numcdf(f, include_weights=True):
    mcdf = make_mcdf(f, mcnubranches)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    if include_weights:
        wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst.geniesyst(f, mcdf.ind)], axis=1) 
        mcdf = pd.concat([mcdf, wgtdf], axis=1)
    return mcdf

def make_mcdf(f, branches=mcbranches, primbranches=mcprimbranches):
    # load the df
    mcdf = loadbranches(f["recTree"], branches)
    
    # Get the correct thing to load
    # Get the correct branches
    while mcdf.columns.nlevels > 2:
        mcdf.columns = mcdf.columns.droplevel(0)

    # No pileup
    # mcdf.index = mcdf.index.droplevel(1)

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
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==211).groupby(level=[0,1]).sum().rename(("npi", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==111).groupby(level=[0,1]).sum().rename(("npi0", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==22).groupby(level=[0,1]).sum().rename(("ng", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==321).groupby(level=[0,1]).sum().rename(("nk", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==310).groupby(level=[0,1]).sum().rename(("nk0", "")))

    # lepton info
    mcdf = mcdf.join(magdf(mcprimdf.genp).groupby(level=[0,1]).first().rename(("lmom", "")))

    # Add in header info
    hdrdf = loadbranches(f["recTree"], hdrbranches)
    hdrdf = hdrdf.rec.hdr

    hproc = hdrdf.proc
    hproc.name = ("proc", "")
    mcdf = mcdf.join(hproc, on=["entry"])
    hcluster = hdrdf.cluster
    hcluster.name = ("cluster", "")
    mcdf = mcdf.join(hcluster, on=["entry"])

    return mcdf

def add_header(f, df):
    def add_levels(df, s):
        return tuple([s] + [""]*(df.columns.nlevels-1))

    hdrdf = loadbranches(f["recTree"], hdrbranches)
    hdrdf = hdrdf.rec.hdr

    evt = hdrdf.evt.copy()
    evt.name = add_levels(df, "evt")
    df = df.join(evt)

    proc = hdrdf.proc.copy()
    proc.name = add_levels(df, "proc")
    df = df.join(proc)

    cluster = hdrdf.cluster.copy()
    cluster.name = add_levels(df, "cluster")
    df = df.join(cluster)

    run = hdrdf.run.copy()
    run.name = add_levels(df, "run")
    df = df.join(run)

    subrun = hdrdf.subrun.copy()
    subrun.name = add_levels(df, "subrun")
    df = df.join(subrun)

    # load the trigger info if it is there
    try:
        triggerdt = hdrdf.triggerinfo.trigger_within_gate.copy()
        triggerdt.name = add_levels(df, "triggerdt")
        df = df.join(triggerdt)
    except:
        pass

    return df

def make_slcdf(f):
    slcdf = loadbranches(f["recTree"], slcbranches)
    slcdf = slcdf.rec

    slc_mcdf = make_mcdf(f, slc_mcbranches, slc_mcprimbranches)
    slc_mcdf.columns = pd.MultiIndex.from_tuples([tuple(["truth"] + list(c)) for c in slc_mcdf.columns])
    slcdf = slcdf.merge(slc_mcdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    slcdf = add_header(f, slcdf)

    trkdf = loadbranches(f["recTree"], trkbranches)
    trkdf = trkdf.rec.slc.reco.pfp.trk[trkdf.rec.slc.reco.pfp.trackScore > 0.5]

    n_proton_reco = ((np.abs(trkdf.truth.p.pdg) == 2212) & (trkdf.truth.p.start_process == 0)).groupby(level=[0,1]).sum()
    n_proton_reco.name = ("truth", "n_proton_reco", "")
    n_proton_reco = n_proton_reco.fillna(0)
    slcdf = slcdf.join(n_proton_reco)

    longest_trk = trkdf.len.groupby(level=[0,1]).max()
    longest_trk.name = ("slc", "longest_trk_len", "")
    longest_trk = longest_trk.fillna(0.)
    slcdf = slcdf.join(longest_trk)

    ntrk = trkdf.len.groupby(level=[0,1]).size()
    ntrk.name = ("slc", "ntrk", "")
    ntrk = ntrk.fillna(0) #.astype(int)
    slcdf = slcdf.join(ntrk)

    return slcdf

def make_stubinfo(f):
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

    def dEdx2dQdx(dEdx, Efield=0.5):
        alpha = 0.930
        rho = 1.38434
        #Efield = 0.5
        beta = 0.212 / (rho * Efield)
        Wion = 1e3 / 4.237e7
        return np.log(alpha + dEdx*beta) / (Wion*beta)
    MIP_dqdx = dEdx2dQdx(1.6) # dEdx2dQdx(1.8, (stubhitdf.efield_vtx + stubhitdf.efield_end) / 2.) * (0.01420 / 1.59e-2)

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

    stub_length_bins = [0, 0.5, 1, 2, 3]
    stub_length_name = ["l0_5cm", "l1cm", "l2cm", "l3cm"]

    dqdx_var = stubdf.inc_sub_charge / stubdf.length
    length_var = stubdf.length

    dqdxs = []
    for blo, bhi, name in zip(stub_length_bins[:-1], stub_length_bins[1:], stub_length_name):
        dqdx = dqdx_var[(length_var > blo) & (length_var < bhi)].groupby(level=[0,1]).max()
        dqdx.name = ("stub", name, "dqdx", "", "") 
        dqdxs.append(dqdx)

    return pd.concat(dqdxs, axis=1)

def make_trkdf(f, include_mc=True):
    # Get trkdf
    trkdf = loadbranches(f["recTree"], trkbranches + trk_mc_branches if include_mc else trkbranches)
    trkdf = trkdf.rec.slc.reco.pfp.trk[trkdf.rec.slc.reco.pfp.trackScore > 0.5]
    trkdf["tindex"] = trkdf.index.get_level_values(2) 

    # Add in truth matching info
    # trkmtchdf = loadbranches(f["recTree"], trkmatchbranches)
    # trkmtchdf = trkmtchdf.rec.slc.reco.pfp.trk.truth.matches
    # trkdf["bestmatch_energy"] = trkmtchdf.energy.groupby(level=[0, 1, 2]).first()

    # Add in daughter track info
    # d = loadbranches(f["recTree"], ["rec.slc.reco.pfp.daughters"]).rec.slc.reco.pfp
    # d.columns = pd.MultiIndex.from_tuples([(b, "", "", "") for b in d.columns])
    # d = d.reset_index().merge(trkdf, right_on=["entry", "rec.slc..index", "id"], 
    #                        left_on=["entry", "rec.slc..index", "daughters"], 
    #                        validate="one_to_one",
    #                        how="inner")
    # longest_d = d.groupby(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"]).len.idxmax()
    # longest_d = d.loc[longest_d].set_index(["entry", "rec.slc..index", "rec.slc.reco.pfp..index"])

    # trkdf[("trkd", "len", "", "")] = longest_d.len
    # trkdf[("trkd", "dir", "x", "")] = longest_d.dir.x
    # trkdf[("trkd", "dir", "y", "")] = longest_d.dir.y
    # trkdf[("trkd", "dir", "z", "")] = longest_d.dir.z

    return trkdf

def make_trunkhitdf(f):
    slcdf = loadbranches(f["recTree"], slcbranches)
    slcdf = slcdf.rec.slc
    trkdf = make_trkdf(f, False)

    trk_distance_to_vtx = dist(slcdf.vertex, trkdf.start).dropna()
    trunk_idx = trk_distance_to_vtx.groupby(level=[0,1]).idxmin()

    return make_trkhitdf(f, trunk_idx)

def make_branchhitdf(f):
    slcdf = loadbranches(f["recTree"], slcbranches)
    slcdf = slcdf.rec.slc
    trkdf = make_trkdf(f, False)

    trk_distance_to_vtx = dist(slcdf.vertex, trkdf.start).dropna()
    trunk_idx = trk_distance_to_vtx.groupby(level=[0,1]).idxmin()
    nottrunk = trkdf.index.difference(trunk_idx).intersection(trk_distance_to_vtx.index)
    branch_idx = trk_distance_to_vtx.loc[nottrunk].groupby(level=[0,1]).idxmin()

    return make_trkhitdf(f, branch_idx)

def make_trkhitdf(f, trkidx):
    trkhitdf = loadbranches(f["recTree"], trkhitbranches)
    trkhitdf = trkhitdf.rec.slc.reco.pfp.trk.calo.I2.points
    idxnames = trkhitdf.index.names

    trkhitdf = trkhitdf.reset_index().set_index(idxnames[:-1])
    trkhitdf = trkhitdf.loc[trkhitdf.index.intersection(trkidx)].reset_index().set_index(idxnames)
    trkhitdf.columns = ["_".join(c).rstrip("_") for c in trkhitdf.columns]
    return add_header(f, trkhitdf)

def make_evtdf(f, verbose=False, load_hits=True, include_mc=True, include_weights=True):
    if verbose: print("starting!")
    slcdf = loadbranches(f["recTree"], slcbranches)
    shwdf = loadbranches(f["recTree"], shwbranches)
    slcdf = slcdf.rec.slc
    shwdf = shwdf.rec.slc.reco.pfp.shw[shwdf.rec.slc.reco.pfp.trackScore < 0.5]
    if verbose: print("df's loaded")

    trkdf = make_trkdf(f, include_mc)
    if verbose: print("made trkdf")

    trk_distance_to_vtx = dist(slcdf.vertex, trkdf.start).dropna()

    if verbose: print("computed distance to vertex")

    # Get the "Trunk" tracks
    trunk = trk_distance_to_vtx.groupby(level=[0,1]).idxmin()
    trunkdf = trkdf.loc[trunk].droplevel(2)

    if verbose: print("load trunk df")

    # Get the "Branch" tracks
    nottrunk = trkdf.index.difference(trunk).intersection(trk_distance_to_vtx.index)

    # branch_dist = magdf((trkdf.start - trunkdf.start)*trunkdf.dir)
    branch = trk_distance_to_vtx.loc[nottrunk].groupby(level=[0,1]).idxmin()
    branchdf = trkdf.loc[branch].droplevel(2)

    if verbose: print("load branch df")

    # Fix column names
    trunkdf.columns = pd.MultiIndex.from_tuples([tuple(["trunk"] + list(c)) for c in trunkdf.columns])
    branchdf.columns = pd.MultiIndex.from_tuples([tuple(["branch"] + list(c)) for c in branchdf.columns])
    evtdf = slcdf.copy()
    evtdf.columns = pd.MultiIndex.from_tuples([tuple(["slc"] + list(c) + ["", ""]) for c in evtdf.columns])

    # Merge it all together
    evtdf = evtdf.merge(trunkdf, left_index=True, right_index=True, how="inner", validate="one_to_one")
    evtdf = evtdf.merge(branchdf, left_index=True, right_index=True, how="inner", validate="one_to_one")

    if verbose: print("merge branch+trunk into evtdf")

    # Add in the hit information at the start of the trunk
    if load_hits:
        trkhitdf = loadbranches(f["recTree"], trkhitbranches)
        trkhitdf = trkhitdf.rec.slc.reco.pfp.trk.calo.I2.points

        if verbose: print("load trk-hit df")

        trunk_min_rr = evtdf.trunk.len - magdf(evtdf.trunk.start - evtdf.branch.start)
        trunk_min_rr.name = ("minrr", "")
        trkhitdf = trkhitdf.join(trunk_min_rr)

        trunk_start_dqdx = trkhitdf[(trkhitdf.rr > trkhitdf.minrr) & (trkhitdf.rr > 5)].groupby(level=[0,1,2]).dqdx.median()
        trunk_start_dqdx.name = ("trunk", "med_start_dqdx", "", "", "")

        evtdf = evtdf.join(trunk_start_dqdx, how="left", 
                           on=["entry", "rec.slc..index", ("trunk", "tindex")])

        if verbose: print("trk-hit vars in evtdf")

    # is there a third track close by?
    not_trunk_branch = nottrunk.difference(branch)
    evtdf["third_trk_dist"] = trk_distance_to_vtx.loc[not_trunk_branch].groupby(level=[0,1]).min()
    # get the third track
    thirdtrk = trk_distance_to_vtx.loc[not_trunk_branch].groupby(level=[0,1]).idxmin()
    thirdtrkdf = trkdf.loc[thirdtrk].droplevel(2)
    evtdf["third_trk_chi2_proton"] = thirdtrkdf.chi2pid.I2.chi2_proton
    evtdf["third_trk_chi2_muon"] = thirdtrkdf.chi2pid.I2.chi2_muon
    evtdf["third_trk_contained"] = TrkInFV(thirdtrkdf.end)

    # What is the longest shower starting close to the vertex?
    shw_distance_to_vtx = dist(slcdf.vertex, shwdf.start)
    evtdf["max_shw_len"] = shwdf.len[shw_distance_to_vtx < 10.].groupby(level=[0,1]).max()

    if verbose: print("extra objects vars in evtdf")

    # Truth matching
    if include_mc:
        slc_mcdf = make_mcdf(f, slc_mcbranches, slc_mcprimbranches)
        slc_mcdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "truth"] + list(c) + [""]) for c in slc_mcdf.columns])
        evtdf = evtdf.merge(slc_mcdf, left_index=True, right_index=True, how="left", validate="one_to_one")

        if include_weights:
            wgtdf = pd.concat([numisyst.numisyst(evtdf.slc.truth.pdg, evtdf.slc.truth.E), geniesyst.geniesyst(f, evtdf.slc.tmatch.idx)], axis=1) 
            wgtdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "wgt"] + list(c) + [""]) for c in wgtdf.columns])
            evtdf = evtdf.merge(wgtdf, left_index=True, right_index=True, how="left", validate="one_to_one")

        if verbose: print("truth matching done")

        # Stubs!
        evtdf = evtdf.join(make_stubinfo(f))

        if verbose: print("stubs done")

    # Sideband-cut
    if SIDEBAND:
        sidebandcut = evtdf.third_trk_contained & (evtdf.third_trk_chi2_proton < 60) & (evtdf.third_trk_chi2_muon > 40)
        evtdf = evtdf[sidebandcut]

    return add_header(f, evtdf)
