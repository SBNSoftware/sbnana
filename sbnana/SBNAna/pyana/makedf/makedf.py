from pyanalib.panda_helpers import *
from .branches import *
from .util import *
from .ttree import *
from . import geniesyst
from . import numisyst


def make_df(f,branches,inds=None):
    tree_name = get_tree_name(f, 'recTree')
    df = loadbranches(f[tree_name], branches)
    if inds:
        df = df.set_index(getcolumns(inds))
        df.index.names = [ind[-1] for ind in getcolumns(inds)] #get bare name of index
        print(df.index.names)
        if '__ntuple' in df.index.names:
            df.reset_index(level='__ntuple', drop=True, inplace=True)
    return df
def make_hdrdf(f,inds=hdrinds):
    df = make_df(f,hdrbranches,inds=inds)
    hdr = df.rec.hdr
    return hdr

def make_mcnudf(f,inds=None, include_weights=False):
    if inds:
        df = make_df(f, mcnubranches+inds,inds=inds)
    else:
        df = make_df(f, mcnubranches)
    mcdf = df.rec.mc.nu
    #----Not supported yet
    if include_weights:
        wgtdf = pd.concat([numisyst.numisyst(mcdf.pdg, mcdf.E), geniesyst.geniesyst(f, mcdf.ind)], axis=1)
        mcdf = pd.concat([mcdf, wgtdf], axis=1)
    return mcdf

def make_trkdf(f, scoreCut=False, requiret0=False, requireCosmic=False):
    tree_name = get_tree_name(f, 'recTree')
    trkdf = loadbranches(f[tree_name], trkbranches)
    if scoreCut:
        trkdf = trkdf.rec.slc.reco[trkdf.rec.slc.reco.pfp.trackScore > 0.5]
    else:
        trkdf = trkdf.rec.slc.reco

    if requiret0:
        trkdf = trkdf[~np.isnan(trkdf.pfp.t0)]

    if requireCosmic:
        trkdf = trkdf[trkdf.pfp.parent == -1]

    trkdf["tindex"] = trkdf.index.get_level_values(2)

    # trk_daughterdf = loadbranches(f[tree_name], pfp_daughter_branch).rec.slc.reco.pfp

    return trkdf

def make_costrkdf(f):
    tree_name = get_tree_name(f, 'recTree')
    trkdf = make_trkdf(f, requiret0=True, requireCosmic=True)
    slcdf = loadbranches(f[tree_name], slcbranches).rec
    return multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

def make_trkhitdf(f):
    tree_name = get_tree_name(f, 'recTree')
    return loadbranches(f[tree_name], trkhitbranches).rec.slc.reco.pfp.trk.calo.I2.points

def make_slc_trkdf(f, trkScoreCut=False, trkDistCut=10., cutClearCosmic=True):
    tree_name = get_tree_name(f, 'recTree')
    # load
    trkdf = make_trkdf(f, trkScoreCut)
    slcdf = loadbranches(f[tree_name], slcbranches).rec

    # merge in tracks
    slcdf = multicol_merge(slcdf, trkdf, left_index=True, right_index=True, how="right", validate="one_to_many")

    # distance from vertex to track start
    slcdf["dist_to_vertex"] = dmagdf(slcdf.slc.vertex, slcdf.pfp.trk.start)

    if trkDistCut > 0:
        slcdf = slcdf[slcdf.dist_to_vertex < trkDistCut]
    if cutClearCosmic:
        slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]

    return slcdf
