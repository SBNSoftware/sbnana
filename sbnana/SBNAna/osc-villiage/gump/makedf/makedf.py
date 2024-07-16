from pyanalib.panda_helpers import *
from .branches import *
from .util import *
from .calo import *
from . import numisyst, g4syst, geniesyst
import uproot

PROTON_MASS = 0.938272
NEUTRON_MASS = 0.939565
MUON_MASS = 0.105658
PION_MASS = 0.139570

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

def make_mcnudf(f, include_weights=False):
    mcdf = make_mcdf(f)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    if include_weights:
        wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst.geniesyst(f, mcdf.ind), g4syst.g4syst(f, mcdf.ind)], axis=1)
        mcdf = multicol_concat(mcdf, wgtdf)
    return mcdf

def make_mchdf(f, include_weights=False):
    mcdf = loadbranches(f["recTree"], mchbranches).rec.mc.prtl
    if include_weights:
        wgtdf = numisyst.numisyst(14, mcdf.E) # TODO: what PDG?
        mcdf = pd.concat([mcdf, wgtdf], axis=1)
    return mcdf

def make_trkdf(f, scoreCut=False, requiret0=False, requireCosmic=False, mcs=False):
    trkdf = loadbranches(f["recTree"], trkbranches + shwbranches)
    if scoreCut:
        trkdf = trkdf.rec.slc.reco[trkdf.rec.slc.reco.pfp.trackScore > 0.5]
    else:
        trkdf = trkdf.rec.slc.reco

    if requiret0:
        trkdf = trkdf[~np.isnan(trkdf.pfp.t0)]

    if requireCosmic:
        trkdf = trkdf[trkdf.pfp.parent == -1]

    if mcs:
        mcsdf = loadbranches(f["recTree"], [trkmcsbranches[0]]).rec.slc.reco.pfp.trk.mcsP
        mcsdf_angle = loadbranches(f["recTree"], [trkmcsbranches[1]]).rec.slc.reco.pfp.trk.mcsP
        mcsdf_angle.index.set_names(mcsdf.index.names, inplace=True)

        mcsdf = mcsdf.merge(mcsdf_angle, how="left", left_index=True, right_index=True)
        mcsgroup = list(range(mcsdf.index.nlevels-1))
        cumlen = mcsdf.seg_length.groupby(level=mcsgroup).cumsum()*14 # convert rad length to cm
        maxlen = (cumlen*(mcsdf.seg_scatter_angles >= 0)).groupby(level=mcsgroup).max()
        trkdf[("pfp", "trk", "mcsP", "len", "", "")] = maxlen


    trkdf[("pfp", "tindex", "", "", "", "")] = trkdf.index.get_level_values(2)

    return trkdf

def make_trkhitdf(f):
    df = loadbranches(f["recTree"], trkhitbranches).rec.slc.reco.pfp.trk.calo.I2.points

    # Firsthit and Lasthit info
    ihit = df.index.get_level_values(-1)
    df["firsthit"] = ihit == 0

    lasthit = df.groupby(level=list(range(df.index.nlevels-1))).tail(1).copy()
    lasthit["lasthit"] = True
    df["lasthit"] = lasthit.lasthit
    df.lasthit = df.lasthit.fillna(False)

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
    while mcdf.columns.nlevels > 2:
        mcdf.columns = mcdf.columns.droplevel(0)

    # Add in primary particle info
    mcprimdf = loadbranches(f["recTree"], primbranches)
    while mcprimdf.columns.nlevels > 2:
        mcprimdf.columns = mcprimdf.columns.droplevel(0)

    mcprimdf.index = mcprimdf.index.rename(mcdf.index.names[:2] + mcprimdf.index.names[2:])

    max_proton_KE = mcprimdf[np.abs(mcprimdf.pdg)==2212].genE.groupby(level=[0,1]).max() - PROTON_MASS
    max_proton_KE.name = ("max_proton_ke", "")
    mcdf = mcdf.join(max_proton_KE)

    mcdf.max_proton_ke = mcdf.max_proton_ke.fillna(0.)

    # particle counts
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

    # particle counts w/ threshold
    proton_KE = mcprimdf[np.abs(mcprimdf.pdg)==2212].genE - PROTON_MASS
    muon_KE = mcprimdf[np.abs(mcprimdf.pdg)==13].genE - MUON_MASS
    pion_KE = mcprimdf[np.abs(mcprimdf.pdg)==211].genE - PION_MASS
    mcdf = mcdf.join(((np.abs(mcprimdf.pdg)==2212) & (proton_KE > 0.05)).groupby(level=[0,1]).sum().rename(("np_50MeV","")))
    mcdf = mcdf.join(((np.abs(mcprimdf.pdg)==2212) & (proton_KE > 0.02)).groupby(level=[0,1]).sum().rename(("np_20MeV","")))
    mcdf = mcdf.join(((np.abs(mcprimdf.pdg)==13) & (muon_KE > 0.02)).groupby(level=[0,1]).sum().rename(("nmu_20MeV","")))
    mcdf = mcdf.join(((np.abs(mcprimdf.pdg)==211) & (pion_KE > 0.04)).groupby(level=[0,1]).sum().rename(("npi_40MeV","")))

 
    # lepton info
    mudf = mcprimdf[np.abs(mcprimdf.pdg)==13].sort_values(mcprimdf.index.names[:2] + [("genE", "")]).groupby(level=[0,1]).last()
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])

    pdf = mcprimdf[mcprimdf.pdg==2212].sort_values(mcprimdf.index.names[:2] + [("genE", "")]).groupby(level=[0,1]).last()
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])

    mcdf = multicol_merge(mcdf, mudf, left_index=True, right_index=True, how="left", validate="one_to_one")
    mcdf = multicol_merge(mcdf, pdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    return mcdf

def make_slc_trkdf(f, trkScoreCut=False, trkDistCut=10., cutClearCosmic=True, **trkArgs):
    # load
    trkdf = make_trkdf(f, trkScoreCut, **trkArgs)
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

def make_eslc_partdf(f, trkDistCut=-1, **trkArgs):
    # load
    partdf = make_epartdf(f, **trkArgs)
    partdf.columns = pd.MultiIndex.from_tuples([tuple(["particle"] + list(c)) for c in partdf.columns])
    eslcdf = make_eslcdf(f)

    # merge in tracks
    eslcdf = multicol_merge(eslcdf, partdf, left_index=True, right_index=True, how="right", validate="one_to_many")
    eslcdf = multicol_add(eslcdf, dmagdf(eslcdf.vertex, eslcdf.particle.start_point).rename("dist_to_vertex"))

    if trkDistCut > 0:
        eslcdf = eslcdf[eslcdf.dist_to_vertex < trkDistCut]

    return eslcdf

def make_stubs(f, det="ICARUS"):
    alpha_sbnd = 0.930                     
    LAr_density_gmL_sbnd = 1.38434
    Efield_sbnd = 0.5                           
    beta_sbnd = 0.212 / (LAr_density_gmL_sbnd * Efield_sbnd)  
    
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
        if det == "SBND":
            return np.log(alpha_sbnd + dEdx*beta_sbnd) / (Wion*beta_sbnd)
        beta = MODB_mc / (LAr_density_gmL_mc * Efield_mc)
        alpha = MODA_mc
        return np.log(alpha + dEdx*beta) / (Wion*beta)
    def dEdx2dQdx_data(dEdx): # data parameters
        if det == "SBND":
            return np.log(alpha_sbnd + dEdx*beta_sbnd) / (Wion*beta_sbnd)
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
    stubdf["truth_pdg"] = stubdf.truth.p.pdg
    stubdf["truth_interaction_id"] = stubdf.truth.p.interaction_id 
    stubdf["truth_gen_E"] = stubdf.truth.p.genE 

    # convert charge to energy
    if ismc:
        stubdf["ke"] = Q2KE_mc(stubdf.Q)
        # also do calorimetric variations
        # TODO: Systematic variations
        stubdf["ke_callo"] = np.nan # Q2KE_mc_callo(stubdf.Q)
        stubdf["ke_calhi"] = np.nan # Q2KE_mc_calhi(stubdf.Q)
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

def make_eslcdf(f):
    eslcdf = loadbranches(f["recTree"], eslcbranches)
    eslcdf = eslcdf.rec.dlp

    etintdf = loadbranches(f["recTree"], etruthintbranches)
    etintdf = etintdf.rec.dlp_true
    
    # match to the truth info
    mcdf = make_mcdf(f)
    # mc is truth
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["truth"] + list(c)) for c in mcdf.columns])

    # Do matching
    # 
    # First get the ML true particle IDs matched to each reco particle
    eslc_matchdf = loadbranches(f["recTree"], eslcmatchedbranches)
    eslc_match_overlap_df = loadbranches(f["recTree"], eslcmatchovrlpbranches)
    eslc_match_overlap_df.index.names = eslc_matchdf.index.names

    eslc_matchdf = multicol_merge(eslc_matchdf, eslc_match_overlap_df, left_index=True, right_index=True, how="left", validate="one_to_one")
    eslc_matchdf = eslc_matchdf.rec.dlp

    # Then use bestmatch.match to get the nu ids in etintdf
    eslc_matchdf_wids = pd.merge(eslc_matchdf, etintdf, left_on=["entry", "match"], right_on=["entry", "id"], how="left")
    eslc_matchdf_wids.index = eslc_matchdf.index

    # Now use nu_ids to get the true interaction information
    eslc_matchdf_trueints = multicol_merge(eslc_matchdf_wids, mcdf, left_on=["entry", "nu_id"], right_index=True, how="left")
    eslc_matchdf_trueints.index = eslc_matchdf_wids.index

    # delete unnecesary matching branches
    del eslc_matchdf_trueints[("match", "")]
    del eslc_matchdf_trueints[("nu_id", "")]
    del eslc_matchdf_trueints[("id", "")]

    # first match is best match
    bestmatch = eslc_matchdf_trueints.groupby(level=list(range(eslc_matchdf_trueints.index.nlevels-1))).first()

    # add extra levels to eslcdf columns
    eslcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) + [""]*2) for c in eslcdf.columns])

    eslcdf_withmc = multicol_merge(eslcdf, bestmatch, left_index=True, right_index=True, how="left")

    # Fix position names (I0, I1, I2) -> (x, y, z)
    def mappos(s):
        if s == "I0": return "x"
        if s == "I1": return "y"
        if s == "I2": return "z"
        return s
    def fixpos(c):
        if c[0] not in ["end_point", "start_point", "start_dir", "vertex", "momentum"]: return c
        return tuple([c[0]] + [mappos(c[1])] + list(c[2:]))

    eslcdf_withmc.columns = pd.MultiIndex.from_tuples([fixpos(c) for c in eslcdf_withmc.columns])

    return eslcdf_withmc

def make_epartdf(f):
    epartdf = loadbranches(f["recTree"], eparticlebranches)
    epartdf = epartdf.rec.dlp.particles

    tpartdf = loadbranches(f["recTree"], trueparticlebranches)
    tpartdf = tpartdf.rec.true_particles
    # cut out EMShowerDaughters
    # tpartdf = tpartdf[(tpartdf.parent == 0)]

    etpartdf = loadbranches(f["recTree"], etrueparticlebranches)
    etpartdf = etpartdf.rec.dlp_true.particles
    etpartdf.columns = [s for s in etpartdf.columns]
    
    # Do matching
    # 
    # First get the ML true particle IDs matched to each reco particle
    epart_matchdf = loadbranches(f["recTree"], eparticlematchedbranches)
    epart_match_overlap_df = loadbranches(f["recTree"], eparticlematchovrlpbranches)
    epart_match_overlap_df.index.names = epart_matchdf.index.names
    epart_matchdf = multicol_merge(epart_matchdf, epart_match_overlap_df, left_index=True, right_index=True, how="left", validate="one_to_one")
    epart_matchdf = epart_matchdf.rec.dlp.particles
    # get the best match (highest match_overlap), assume it's sorted
    bestmatch = epart_matchdf.groupby(level=list(range(epart_matchdf.index.nlevels-1))).first()
    bestmatch.columns = [s for s in bestmatch.columns]

    # Then use betmatch.match to get the G4 track IDs in etpartdf
    bestmatch_wids = pd.merge(bestmatch, etpartdf, left_on=["entry", "match"], right_on=["entry", "id"], how="left")
    bestmatch_wids.index = bestmatch.index

    # Now use the G4 track IDs to get the true particle information
    bestmatch_trueparticles = multicol_merge(bestmatch_wids, tpartdf, left_on=["entry", "track_id"], right_on=["entry", ("G4ID", "")], how="left")
    bestmatch_trueparticles.index = bestmatch_wids.index

    # delete unnecesary matching branches
    del bestmatch_trueparticles[("match", "")]
    del bestmatch_trueparticles[("track_id", "")]
    del bestmatch_trueparticles[("id", "")]

    # add extra level to epartdf columns
    epartdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) + [""]) for c in epartdf.columns])

    # put everything in epartdf
    for c in bestmatch_trueparticles.columns:
        epartdf[tuple(["truth"] + list(c))] = bestmatch_trueparticles[c]

    # Fix position names (I0, I1, I2) -> (x, y, z)
    def mappos(s):
        if s == "I0": return "x"
        if s == "I1": return "y"
        if s == "I2": return "z"
        return s
    def fixpos(c):
        if c[0] not in ["end_point", "start_point", "start_dir", "vertex"]: return c
        return tuple([c[0]] + [mappos(c[1])] + list(c[2:]))

    epartdf.columns = pd.MultiIndex.from_tuples([fixpos(c) for c in epartdf.columns])

    return epartdf

def make_eevtdf(f):
    # load slices and particles
    partdf = make_epartdf(f)

    df = make_eslcdf(f)

    # load the proton and muon candidates
    primary = partdf.is_primary
    mudf = partdf[primary & (partdf.pid == 2)].sort_values(partdf.index.names[:2] + [("length", "", "")]).groupby(level=[0,1]).last()
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])

    pdf = partdf[primary & (partdf.pid == 4)].sort_values(partdf.index.names[:2] + [("length", "", "")]).groupby(level=[0,1]).last()
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])

    df = multicol_merge(df, mudf, left_index=True, right_index=True, how="left", validate="one_to_one")
    df = multicol_merge(df, pdf, left_index=True, right_index=True, how="left", validate="one_to_one")

    # in case we want to cut out other objects -- save the highest energy of each other particle
    lead_gamma_energy = partdf.ke[primary & (partdf.pid == 0)].groupby(level=[0,1]).max().rename("lead_gamma_energy")
    df = multicol_add(df, lead_gamma_energy)

    lead_elec_energy = partdf.ke[primary & (partdf.pid == 1)].groupby(level=[0,1]).max().rename("lead_elec_energy")
    df = multicol_add(df, lead_elec_energy)

    lead_pion_length = partdf.length[primary & (partdf.pid == 3)].groupby(level=[0,1]).max().rename("lead_pion_length")
    df = multicol_add(df, lead_pion_length)

    subl_muon_length = partdf[primary & (partdf.pid == 2)].sort_values(partdf.index.names[:2] + [("length", "", "")]).length.groupby(level=[0,1]).nth(-2).rename("subl_muon_length")
    df = multicol_add(df, subl_muon_length)

    subl_proton_length = partdf[primary & (partdf.pid == 4)].sort_values(partdf.index.names[:2] + [("length", "", "")]).length.groupby(level=[0,1]).nth(-2).rename("subl_proton_length")
    df = multicol_add(df, subl_proton_length)

    # Apply pre-selection: Require fiducial vertex, at least one muon, at least one proton

    # require both muon and proton to be present
    df = df[~np.isnan(df.mu.pid) & ~np.isnan(df.p.pid)]

    # require fiducial verex
    df = df[InFV(df.vertex, 50)]

    return df


def make_evtdf(f, trkScoreCut=False, trkDistCut=10., cutClearCosmic=True, **trkArgs):
    
    # ----- sbnd or icarus? -----
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"
    
    mcdf = make_mcdf(f)
    trkdf = make_trkdf(f, trkScoreCut, **trkArgs)
    slcdf = make_slcdf(f)
    stubdf = make_stubs(f, det=DETECTOR)
    
    # ----- merge dfs -----
    # load stubs
    slcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +[""]) for c in slcdf.columns])
    slcdf = slcdf.join(stubdf)
    # load pfps
    # slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    trkdf = multicol_add(trkdf, dmagdf(slcdf.slc.vertex, trkdf.pfp.trk.start).rename(("pfp", "dist_to_vertex")))
    if trkDistCut > 0:
        trkdf = trkdf[trkdf.pfp.dist_to_vertex < trkDistCut]
    if cutClearCosmic:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

    # ---- calculate additional info ----
    
    # track containment
    trkdf[("pfp", "trk", "is_contained", "", "", "")] = (InFV(trkdf.pfp.trk.start, 0, det=DETECTOR)) & (InFV(trkdf.pfp.trk.end, 0, det=DETECTOR))

    # reco momentum -- range for contained, MCS for exiting
    trkdf[("pfp", "trk", "P", "p_muon", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_muon", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_muon", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_muon","", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_muon", "", "")]

    trkdf[("pfp", "trk", "P", "p_pion", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_pion", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_pion", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_pion", "", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_pion", "", "")]

    trkdf[("pfp", "trk", "P", "p_proton", "", "")] = np.nan
    trkdf.loc[trkdf.pfp.trk.is_contained, ("pfp", "trk", "P", "p_proton", "", "")]  = trkdf.loc[(trkdf.pfp.trk.is_contained), ("pfp", "trk", "rangeP", "p_proton", "", "")]
    trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "P", "p_proton", "", "")] = trkdf.loc[np.invert(trkdf.pfp.trk.is_contained), ("pfp", "trk", "mcsP", "fwdP_proton", "", "")]

    # opening angles
    trkdf[("pfp", "trk", "cos", "x", "", "")] = np.nan
    trkdf[("pfp", "trk", "cos", "x", "", "")] = (trkdf.pfp.trk.end.x-trkdf.pfp.trk.start.x)/trkdf.pfp.trk.len
    trkdf[("pfp", "trk", "cos", "y", "", "")] = np.nan
    trkdf[("pfp", "trk", "cos", "y", "", "")] = (trkdf.pfp.trk.end.y-trkdf.pfp.trk.start.y)/trkdf.pfp.trk.len
    trkdf[("pfp", "trk", "cos", "z", "", "")] = np.nan
    trkdf[("pfp", "trk", "cos", "z", "", "")] = (trkdf.pfp.trk.end.z-trkdf.pfp.trk.start.z)/trkdf.pfp.trk.len

    # ----- loose PID for candidates ----
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = np.nan
    trkdf[("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")] = trkdf.pfp.trk.chi2pid.I2.chi2_muon/trkdf.pfp.trk.chi2pid.I2.chi2_proton
    
    # mu candidate is track pfp with smallest chi2_mu/chi2_p
    mudf = trkdf[(trkdf.pfp.trackScore > 0.5)].sort_values(trkdf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).head(1)
    mudf.columns = pd.MultiIndex.from_tuples([tuple(["mu"] + list(c)) for c in mudf.columns])
    slcdf = multicol_merge(slcdf, mudf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_mu = mudf.index
        
    # p candidate is track pfp with largest chi2_mu/chi2_p of remaining pfps
    idx_pfps = trkdf.index
    idx_not_mu = idx_pfps.difference(idx_mu)
    notmudf = trkdf.loc[idx_not_mu]
    pdf = notmudf[(notmudf.pfp.trackScore > 0.5)].sort_values(notmudf.pfp.index.names[:-1] + [("pfp", "trk", "chi2pid", "I2", "mu_over_p", "")]).groupby(level=[0,1]).tail(1)
    pdf.columns = pd.MultiIndex.from_tuples([tuple(["p"] + list(c)) for c in pdf.columns])
    slcdf = multicol_merge(slcdf, pdf.droplevel(-1), left_index=True, right_index=True, how="left", validate="one_to_one")
    idx_p = pdf.index

    # calculated and save transverse kinematic variables
    mudf = slcdf.mu.pfp.trk
    pdf = slcdf.p.pfp.trk

    # Caculate transverse kinematics
    MASS_A = 22*NEUTRON_MASS + 18*PROTON_MASS - 0.34381
    BE = 0.0295
    MASS_Ap = MASS_A - NEUTRON_MASS + BE

    mu_p = mudf.P.p_muon
    mu_p_x = mu_p * mudf.cos.x
    mu_p_y = mu_p * mudf.cos.y
    mu_p_z = mu_p * mudf.cos.z
    mu_phi_x = mu_p_x/mag2d(mu_p_x, mu_p_y)
    mu_phi_y = mu_p_y/mag2d(mu_p_x, mu_p_y)

    p_p = pdf.P.p_proton
    p_p_x = p_p * pdf.cos.x
    p_p_y = p_p * pdf.cos.y
    p_p_z = p_p * pdf.cos.z
    p_phi_x = p_p_x/mag2d(p_p_x, p_p_y)
    p_phi_y = p_p_y/mag2d(p_p_x, p_p_y)

    mu_Tp_x = mu_phi_y*mu_p_x - mu_phi_x*mu_p_y
    mu_Tp_y = mu_phi_x*mu_p_x - mu_phi_y*mu_p_y
    mu_Tp = mag2d(mu_Tp_x, mu_Tp_y)

    p_Tp_x = mu_phi_y*p_p_x - mu_phi_x*p_p_y
    p_Tp_y = mu_phi_x*p_p_x - mu_phi_y*p_p_y
    p_Tp = mag2d(p_Tp_x, p_Tp_y)

    del_Tp_x = mu_Tp_x + p_Tp_x
    del_Tp_y = mu_Tp_y + p_Tp_y
    del_Tp = mag2d(del_Tp_x, del_Tp_y)

    del_alpha = np.arccos(-(mu_Tp_x*del_Tp_x + mu_Tp_y*del_Tp_y)/(mu_Tp*del_Tp))
    del_phi = np.arccos(-(mu_Tp_x*p_Tp_x + mu_Tp_y*p_Tp_y)/(mu_Tp*p_Tp))

    mu_E = mag2d(mu_p, MUON_MASS)
    p_E = mag2d(p_p, PROTON_MASS)

    R = MASS_A + mu_p_z + p_p_z - mu_E - p_E
    del_Lp = 0.5*R - mag2d(MASS_Ap, del_Tp)**2/(2*R)
    del_p = mag2d(del_Tp, del_Lp)

    del_alpha = del_alpha.rename("del_alpha")
    slcdf = multicol_add(slcdf, del_alpha)
    del_phi = del_phi.rename("del_phi")
    slcdf = multicol_add(slcdf, del_phi)
    del_Tp = del_Tp.rename("del_Tp")
    slcdf = multicol_add(slcdf, del_Tp)
    del_p = del_p.rename("del_p")
    slcdf = multicol_add(slcdf, del_p)
    
    # note if there are any other track/showers
    idx_not_mu_p = idx_not_mu.difference(idx_p)
    otherdf = trkdf.loc[idx_not_mu_p]
    # longest other shower
    othershwdf = otherdf[otherdf.pfp.trackScore < 0.5]
    other_shw_length = othershwdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_shw_length")
    slcdf = multicol_add(slcdf, other_shw_length)
    # longest other track
    othertrkdf = otherdf[otherdf.pfp.trackScore > 0.5]
    other_trk_length = othertrkdf.pfp.trk.len.groupby(level=[0,1]).max().rename("other_trk_length")
    slcdf = multicol_add(slcdf, other_trk_length)
 
    # ---- truth match ----
    bad_tmatch = np.invert(slcdf.slc.tmatch.eff > 0.5) & (slcdf.slc.tmatch.idx >= 0)
    slcdf.loc[bad_tmatch, ("slc","tmatch","idx", "", "", "", "")] = np.nan

    mcdf.columns = pd.MultiIndex.from_tuples([tuple(list(c) +["", "", "", "", ""]) for c in mcdf.columns])     # match # of column levels

    # require both muon and proton to be present
    slcdf = slcdf[~np.isnan(mudf.P.p_muon) & ~np.isnan(pdf.P.p_muon)]

    df = multicol_merge(slcdf.reset_index(), 
                  mcdf.reset_index(),
                  left_on=[("entry", "", "",), 
                           ("slc", "tmatch", "idx")], 
                  right_on=[("entry", "", ""), 
                            ("rec.mc.nu..index", "", "")], 
                  how="left"
                  ) 

    df = df.set_index(slcdf.index.names, verify_integrity=True)
    
    return df
