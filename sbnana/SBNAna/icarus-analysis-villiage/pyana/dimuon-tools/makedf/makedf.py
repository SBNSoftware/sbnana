from pyanalib.panda_helpers import *
from .branches import *
from .util import *
from . import geniesyst
from . import numisyst

def make_hdrdf(f):
    hdr = loadbranches(f["recTree;1"], hdrbranches).rec.hdr
    return hdr

def make_mcnudf(f, include_weights=True):
    mcdf = make_mcdf(f)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    if include_weights:
        wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst.geniesyst(f, mcdf.ind)], axis=1)
        mcdf = pd.concat([mcdf, wgtdf], axis=1)
    return mcdf

def make_trkdf(f, scoreCut=False, requiret0=False, requireCosmic=False):
    trkdf = loadbranches(f["recTree"], trkbranches + shwbranches)
    if scoreCut:
        trkdf = trkdf.rec.slc.reco[trkdf.rec.slc.reco.pfp.trackScore > 0.5]
    else:
        trkdf = trkdf.rec.slc.reco

    if requiret0:
        trkdf = trkdf[~np.isnan(trkdf.pfp.t0)]

    if requireCosmic:
        trkdf = trkdf[trkdf.pfp.parent == -1]

    trkdf[("pfp", "tindex", "", "", "", "")] = trkdf.index.get_level_values(2)

    # trk_daughterdf = loadbranches(f["recTree;1"], pfp_daughter_branch).rec.slc.reco.pfp

    return trkdf

def make_costrkdf(f):
    trkdf = make_trkdf(f, requiret0=True, requireCosmic=True)
    slcdf = loadbranches(f["recTree;1"], slcbranches).rec
    return multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

def make_trkhitadcdf(f):
    return loadbranches(f["recTree"], trkhitadcbranches).rec.slc.reco.pfp.trk.calo.I2.points

def make_trkhitdf(f, include_adc=False):
    df = loadbranches(f["recTree"], trkhitbranches).rec.slc.reco.pfp.trk.calo.I2.points

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
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==211).groupby(level=[0,1]).sum().rename(("npi", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==111).groupby(level=[0,1]).sum().rename(("npi0", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==22).groupby(level=[0,1]).sum().rename(("ng", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==321).groupby(level=[0,1]).sum().rename(("nk", "")))
    mcdf = mcdf.join((np.abs(mcprimdf.pdg)==310).groupby(level=[0,1]).sum().rename(("nk0", "")))

    # lepton info
    mcdf = mcdf.join(magdf(mcprimdf.genp).groupby(level=[0,1]).first().rename(("lmom", "")))

    return mcdf

def make_slc_trkdf(f, trkScoreCut=False, trkDistCut=10., cutClearCosmic=True):
    # load
    trkdf = make_trkdf(f, trkScoreCut)
    slcdf = make_slcdf(f)

    # merge in tracks
    slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    # distance from vertex to track start
    slcdf["dist_to_vertex"] = dmagdf(slcdf.slc.vertex, slcdf.pfp.trk.start)

    if trkDistCut > 0:
        slcdf = slcdf[slcdf.dist_to_vertex < trkDistCut]
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

    trunk = evtdf.dist_to_vertex[is_muon(evtdf)].groupby(level=list(range(evtdf.index.nlevels-1))).idxmin()
    trunkdf = evtdf.loc[trunk].droplevel(evtdf.index.nlevels-1)

    nottrunk = evtdf.index.difference(trunk)
    branch = evtdf.dist_to_vertex.loc[nottrunk][is_muon(evtdf)].groupby(level=list(range(evtdf.index.nlevels-1))).idxmin()
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

    # only take collection plane
    stubdf = stubdf[stubdf.plane == 2]

    stub_length_bins = [0, 0.5, 1, 2, 3]
    stub_length_name = ["l0_5cm", "l1cm", "l2cm", "l3cm"]

    dqdx_var = stubdf.inc_sub_charge / stubdf.length
    length_var = stubdf.length

    dqdxs = []
    for blo, bhi, name in zip(stub_length_bins[:-1], stub_length_bins[1:], stub_length_name):
        dqdx = dqdx_var[(length_var > blo) & (length_var < bhi)].groupby(level=[0,1]).max()
        dqdx.name = ("stub", name, "dqdx", "", "", "")
        dqdxs.append(dqdx)

    return pd.concat(dqdxs, axis=1)

def make_evtdf(f, load_hits=True):
    # implement the sideband selection
    def is_muon_contained(trk):
        return (trk.chi2pid.I2.chi2_proton > 30) & (trk.chi2pid.I2.chi2_muon < 80) & TrkInFV(trk.end)
    def is_muon_exiting(trk):
        return (trk.len > 100) & ~TrkInFV(trk.end)
    def is_trk(pfp):
        return pfp.trackScore > 0.5
    def is_muon(pfp):
        return is_trk(pfp) & (is_muon_contained(pfp.trk) | is_muon_exiting(pfp.trk))
    def evt_presel(df):
        trklevel = list(range(df.index.nlevels-1))
        nmuon = broadcast(is_muon(df.pfp).groupby(level=trklevel).sum(), df)
        return SlcInFV(df.slc.vertex) & (nmuon >= 2)

    trk_slcdf = make_slc_trkdf(f)
    trk_slcdf = trk_slcdf[evt_presel(trk_slcdf)]
    slcdf = trk_slcdf.copy()

    trklevel = list(range(slcdf.index.nlevels-1))

    trunk = slcdf.dist_to_vertex[is_muon(slcdf.pfp)].groupby(level=trklevel).idxmin()
    trunkdf = slcdf.loc[trunk].droplevel(slcdf.index.nlevels-1)

    nottrunk = slcdf.index.difference(trunk)
    branch = slcdf.dist_to_vertex.loc[nottrunk][is_muon(slcdf.pfp)].groupby(level=trklevel).idxmin()
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
        trunk_min_rr.name = "minrr" #("minrr", "")
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
    evtdf["third_trk_dist"] = trk_slcdf.dist_to_vertex.loc[not_trunk_branch].groupby(level=trklevel).min()
    # get the third track
    thirdtrk = trk_slcdf.dist_to_vertex.loc[not_trunk_branch].groupby(level=trklevel).idxmin()
    thirdtrkdf = trk_slcdf.loc[thirdtrk].droplevel(len(trklevel)).pfp.trk
    evtdf["third_trk_chi2_proton"] = thirdtrkdf.chi2pid.I2.chi2_proton
    evtdf["third_trk_chi2_muon"] = thirdtrkdf.chi2pid.I2.chi2_muon
    evtdf["third_trk_contained"] = TrkInFV(thirdtrkdf.end)

    # Get the other pfp df
    othrdf = trk_slcdf.loc[not_trunk_branch]
    shwdf = othrdf[othrdf.pfp.trackScore < 0.5]

    # min/max chi2s
    evtdf["min_othr_chi2_proton"] = othrdf.pfp.trk.chi2pid.I2.chi2_proton.groupby(level=trklevel).min()
    evtdf["max_othr_chi2_muon"] = othrdf.pfp.trk.chi2pid.I2.chi2_muon.groupby(level=trklevel).max()
    evtdf["max_shw_len"] = shwdf.pfp.shw.len.groupby(level=trklevel).max()

    # stubs!
    evtdf = evtdf.join(make_stubs(f))

    return evtdf

    
    
