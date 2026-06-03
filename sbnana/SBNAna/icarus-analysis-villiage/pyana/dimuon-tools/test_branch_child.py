from pyanalib.ntuple_glob import NTupleGlob
from pyanalib.panda_helpers import *

fname = "/pnfs/icarus/scratch/users/gputnam/DMCP2023G/jamie-gen/MCNuTemplateD/73445101_999/prodcorsika_standard_icarus_20250122T004832-AddCorsika_20250122T020643-G4_20250122T053515-DetSim_20250122T065127-MCstage0_20250122T073523-MCstage1.flat.caf.root"

ntuples = NTupleGlob(fname, None)

pfpbranch = "rec.slc.reco.pfp."
trkbranch = pfpbranch + "trk."
shwbranch = pfpbranch + "shw."

pfp_daughter_branch = [
    pfpbranch + "daughters"
]

trkbranches = [
    pfpbranch + "id",
    pfpbranch + "ndaughters",
    pfpbranch + "trackScore",
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
]

def f(uf):
   # load dataframes
   trkdf = loadbranches(uf["recTree"], trkbranches).rec.slc.reco.pfp
   ddf  = loadbranches(uf["recTree"], pfp_daughter_branch).rec.slc.reco.pfp

   # get columns to merge with. If you want to save more stuff, put it here
   trkmerge = trkdf[[("id", "", "", ""), 
                     ("trk", "len", "", ""), 
                     ("trackScore", "", "", ""),
                     ("trk", "start", "x", ""),
                     ("trk", "start", "y", ""),
                     ("trk", "start", "z", "")
                    ]]
   # fix column names
   trkmerge.columns = ["id", 
                       "len", 
                       "trackScore",
                       "start_x",
                       "start_y",
                       "start_z"
                      ]

   # Get the length of each daughter track. Maintain the index
   ddf_index = ddf.index.names
   ddf = ddf.reset_index().merge(trkmerge.reset_index(), how="left", suffixes=("", "_t"),
       left_on=["entry", "rec.slc..index", "daughters"], 
      right_on=["entry", "rec.slc..index", "id"]).set_index(ddf_index)

   # select the longest daughter
   # Get the ID of the longest daughter
   longest_idx = ddf.len.groupby(level=[0,1,2]).idxmax()
   longest_ddf = ddf.loc[longest_idx].droplevel(-1)

   # Put information back into the trkdf

   # Put the "longest daughter track length" into trkdf
   # You could also now put in other variables from the longest daughter. Just put them into "trkmerge" above first.
   # Note: in the case of no daughters, these varariables will be NaN
   trkdf["longest_daughter_length"] = longest_ddf.len
   trkdf["longest_daughter_score"] = longest_ddf.trackScore
   trkdf["longest_daughter_start_x"] = longest_ddf.start_x
   trkdf["longest_daughter_start_y"] = longest_ddf.start_y
   trkdf["longest_daughter_start_z"] = longest_ddf.start_z
   return trkdf

dfs = ntuples.dataframes(nproc="auto", fs=[f])
