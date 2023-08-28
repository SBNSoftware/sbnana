hdrbranches = [
    "rec.hdr.pot",
    "rec.hdr.first_in_subrun",
    "rec.hdr.run",
    "rec.hdr.subrun",
    "rec.hdr.ngenevt",
    "rec.hdr.evt",
    "rec.hdr.proc",
    "rec.hdr.cluster",
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

pfpbranch = "rec.slc.reco.pfp."
trkbranch = pfpbranch + "trk."
shwbranch = pfpbranch + "shw."

pfobranches = [
    pfpbranch + "pfochar.chgendfrac",
    pfpbranch + "pfochar.chgfracspread",
    pfpbranch + "pfochar.linfitdiff",
    pfpbranch + "pfochar.linfitlen",
    pfpbranch + "pfochar.linfitgaplen",
    pfpbranch + "pfochar.linfitrms",
    pfpbranch + "pfochar.openanglediff",
    pfpbranch + "pfochar.pca2ratio",
    pfpbranch + "pfochar.pca3ratio", 
    pfpbranch + "pfochar.vtxdist" 
]

pfpbranches = [
    pfpbranch + "parent_is_primary",
    pfpbranch + "slcID",
    pfpbranch + "trackScore",
    pfpbranch + "parent",
    pfpbranch + "id",
    pfpbranch + "t0",
] + pfobranches

pfp_daughter_branch = [
    pfpbranch + "daughters"
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
    trkbranch + "chi2pid.2.pid_ndof",
    trkbranch + "chi2pid.2.chi2_muon",
    trkbranch + "chi2pid.2.chi2_proton",
    trkbranch + "chi2pid.2.pida",
] + pfpbranches

shwbranches = [
  shwbranch + "len"
]

trkhitadcbranches = [
  trkbranch + "calo.2.points.adcs"
]

trkhitbranches = [
    trkbranch + "calo.2.points.dedx",
    trkbranch + "calo.2.points.dqdx",
    trkbranch + "calo.2.points.pitch",
    trkbranch + "calo.2.points.integral",
    trkbranch + "calo.2.points.rr",
    trkbranch + "calo.2.points.wire",
    trkbranch + "calo.2.points.tpc",
    trkbranch + "calo.2.points.sumadc",
    trkbranch + "calo.2.points.t",
    trkbranch + "calo.2.points.x",
    trkbranch + "calo.2.points.y",
    trkbranch + "calo.2.points.z",

    #trkbranch + "calo.2.points.width",
    #trkbranch + "calo.2.points.mult",
    #trkbranch + "calo.2.points.tdc0",

    #trkbranch + "calo.2.points.truth.h_e",
    #trkbranch + "calo.2.points.truth.h_nelec",
    #trkbranch + "calo.2.points.truth.pitch",
    #trkbranch + "calo.2.points.truth.rr",
]

for n in trueparticlenames: trkbranches.append(trkbranch + "truth.p." + n)

slcbranches = [
    "rec.slc.is_clear_cosmic",
    "rec.slc.vertex.x", "rec.slc.vertex.y", "rec.slc.vertex.z",
    "rec.slc.self",
    "rec.slc.tmatch.eff",
    "rec.slc.tmatch.pur",
    "rec.slc.tmatch.index",
    "rec.slc.producer",
    "rec.slc.nuid.crlongtrkdiry"
]

mcbranches = [
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

slc_mcbranches = ["rec.slc.truth." + ".".join(s.split(".")[3:]) for s in mcbranches]
slc_mcprimbranches = ["rec.slc.truth." + ".".join(s.split(".")[3:]) for s in mcprimbranches]

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


