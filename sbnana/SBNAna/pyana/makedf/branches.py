slcbranch = "rec.slc."
pfpbranch = slcbranch + "reco.pfp."
trkbranch = pfpbranch + "trk."
shwbranch = pfpbranch + "shw."
hdrbranch = 'rec.hdr.'

hdrinds = [
    hdrbranch + 'run',
    hdrbranch + 'subrun',
    hdrbranch + 'evt',
]

hdrbranches = [
    hdrbranch + "pot",
    hdrbranch + "first_in_subrun",
    hdrbranch + "run",
    hdrbranch + "subrun",
    hdrbranch + "ngenevt",
    hdrbranch + "evt",
    hdrbranch + "proc",
    hdrbranch + "cluster",
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
]

shwbranches = [
 shwbranch + 'bestplane',
 shwbranch + 'bestplane_dEdx',
 shwbranch + 'bestplane_energy',
 shwbranch + 'conversion_gap',
 shwbranch + 'density',
 shwbranch + 'dir.x',
 shwbranch + 'dir.y',
 shwbranch + 'dir.z',
 shwbranch + 'end.x',
 shwbranch + 'end.y',
 shwbranch + 'end.z',
 shwbranch + 'len',
 shwbranch + 'open_angle',
 shwbranch + 'plane.0.dEdx',
 shwbranch + 'plane.0.energy',
 shwbranch + 'plane.0.nHits',
 shwbranch + 'plane.0.wirePitch',
 shwbranch + 'plane.1.dEdx',
 shwbranch + 'plane.1.energy',
 shwbranch + 'plane.1.nHits',
 shwbranch + 'plane.1.wirePitch',
 shwbranch + 'plane.2.dEdx',
 shwbranch + 'plane.2.energy',
 shwbranch + 'plane.2.nHits',
 shwbranch + 'plane.2.wirePitch',
 shwbranch + 'razzle.bestScore',
 shwbranch + 'razzle.electronScore',
 shwbranch + 'razzle.otherScore',
 shwbranch + 'razzle.pdg',
 shwbranch + 'razzle.photonScore',
 shwbranch + 'start.x',
 shwbranch + 'start.y',
 shwbranch + 'start.z',
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
  trkbranch + "chi2pid.2.chi2_muon",
  trkbranch + "chi2pid.2.chi2_proton",
  trkbranch + "chi2pid.2.chi2_kaon",
  trkbranch + "chi2pid.2.chi2_pion",
  trkbranch + 'crttrack.angle',
  trkbranch + 'crttrack.time',
  trkbranch + 'dazzle.bestScore',
  trkbranch + 'dazzle.muonScore',
  trkbranch + 'dazzle.otherScore',
  trkbranch + 'dazzle.pdg',
  trkbranch + 'dazzle.pionScore',
  trkbranch + 'dazzle.protonScore',
]

for n in trueparticlenames: trkbranches.append(trkbranch + "truth.p." + n)
for n in trueparticlenames: shwbranches.append(shwbranch + "truth.p." + n)

pfpallbranches = pfpbranches+trkbranches+shwbranches
trkbranches += pfpbranches
shwbranches += pfpbranches

#This causes a length mismatch
# pfpinds = hdrinds + [
#     pfpbranch + 'slcID',
#     pfpbranch + 'id'
#     ]

pfp_daughter_branch = [
    pfpbranch + "daughters"
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

mcnubranches = [
    "rec.mc.nu.E",
    "rec.mc.nu.position.x",
    "rec.mc.nu.position.y",
    "rec.mc.nu.position.z",
    "rec.mc.nu.pdg",
    "rec.mc.nu.iscc",
    "rec.mc.nu.genie_mode"
]


