from pyanalib.panda_helpers import *
from .branches import *
from .util import *
from . import geniesyst
from . import numisyst

def make_hdrdf(f):
    hdr = loadbranches(f["recTree;1"], hdrbranches).rec.hdr
    return hdr

def make_mcnudf(f, include_weights=True):
    mcdf = loadbranches(f["recTree;1"], mcnubranches).rec.mc.nu
    mcdf["ind"] = mcdf.index.get_level_values(1)
    if include_weights:
        wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst.geniesyst(f, mcdf.ind)], axis=1)
        mcdf = pd.concat([mcdf, wgtdf], axis=1)
    return mcdf

def make_trkdf(f, scoreCut=False, requiret0=False, requireCosmic=False):
    print(all(trkbranches in  f['recTree;1'].keys()))
    trkdf = loadbranches(f["recTree;1"], trkbranches)
    if scoreCut:
        trkdf = trkdf.rec.slc.reco[trkdf.rec.slc.reco.pfp.trackScore > 0.5]
    else:
        trkdf = trkdf.rec.slc.reco

    if requiret0:
        trkdf = trkdf[~np.isnan(trkdf.pfp.t0)]

    if requireCosmic:
        trkdf = trkdf[trkdf.pfp.parent == -1]

    trkdf["tindex"] = trkdf.index.get_level_values(2)

    # trk_daughterdf = loadbranches(f["recTree;1"], pfp_daughter_branch).rec.slc.reco.pfp

    return trkdf

def make_costrkdf(f):
    trkdf = make_trkdf(f, requiret0=True, requireCosmic=True)
    slcdf = loadbranches(f["recTree;1"], slcbranches).rec
    return multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

def make_trkhitdf(f):
    return loadbranches(f["recTree;1"], trkhitbranches).rec.slc.reco.pfp.trk.calo.I2.points

def make_slc_trkdf(f, trkScoreCut=False, trkDistCut=10., cutClearCosmic=True):
    # load
    trkdf = make_trkdf(f, trkScoreCut)
    slcdf = loadbranches(f["recTree;1"], slcbranches).rec

    # merge in tracks
    slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    # distance from vertex to track start
    slcdf["dist_to_vertex"] = dmagdf(slcdf.slc.vertex, slcdf.pfp.trk.start)

    if trkDistCut > 0:
        slcdf = slcdf[slcdf.dist_to_vertex < trkDistCut]
    if cutClearCosmic:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

    return slcdf
