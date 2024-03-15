import numpy as np
import pandas as pd

# CATHODE CROSSING
correction_file = "/exp/icarus/app/users/jzettle/BSM_Systs/Weights_Output/dir_x_weights_cathode.txt"
cathode_corrections = np.loadtxt(correction_file).T
cathode_correction_W_split = cathode_corrections[1]
cathode_correction_W_reco = cathode_corrections[3]
cathode_correction_E_split = cathode_corrections[9]
cathode_correction_E_reco = cathode_corrections[11]
cathode_correction_bins = np.linspace(-1, 1, cathode_corrections.shape[1]+1)
cathode_correction_err_W = 0.3
cathode_correction_err_E = 0.2
def tpc_index(x):
    itpc = (x*0).fillna(0).astype(int)
    itpc[x > -210.29] += 1
    itpc[x > 61.94] += 1
    itpc[x > 210.29] += 1
    return itpc

# Z=0 INDUCTION GAP
correction_file = "/exp/icarus/app/users/jzettle/BSM_Systs/Weights_Output/dir_z_weights_across_z.txt"
gap_corrections = np.loadtxt(correction_file).T
gap_correction_W_split = gap_corrections[1]
gap_correction_W_reco = gap_corrections[2]
gap_correction_E_split = gap_corrections[3]
gap_correction_E_reco = gap_corrections[4]
# don't actually apply corrections for split -- we can't fix these
gap_correction_W_split[:] = 1
gap_correction_E_split[:] = 1
# large error since this is harder to trust
gap_correction_err = 1
gap_correction_bins = np.linspace(-1, 1, cathode_corrections.shape[1]+1)

def correct(evtdf):
    # CATHODE CROSSING
    cathode_correction_trunk_var = evtdf.trunk.trk.dir.x
    cathode_correction_branch_var = evtdf.branch.trk.dir.x

    # TRUNK
    true_trunk_crosses_cathode = tpc_index(evtdf.trunk.trk.truth.p.start.x) != tpc_index(evtdf.trunk.trk.truth.p.end.x)
    true_same_cryo = (tpc_index(evtdf.trunk.trk.truth.p.start.x)//2) == (tpc_index(evtdf.trunk.trk.truth.p.end.x)//2)
    true_trunk_crosses_cathode = true_trunk_crosses_cathode & true_same_cryo
    reco_trunk_crosses_cathode = tpc_index(evtdf.trunk.trk.start.x) != tpc_index(evtdf.trunk.trk.end.x)
    cryoE = evtdf.trunk.trk.producer == 0
    
    # CORRECTION
    cathode_cross_trunk_weight = pd.Series(1, evtdf.index)
    cathode_cross_trunk_weight[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & ~cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & ~cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_W_split)\
          .astype(float).fillna(0)
    cathode_cross_trunk_weight[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & ~cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & ~cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_W_reco)\
          .astype(float).fillna(0)
    cathode_cross_trunk_weight[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_E_split)\
          .astype(float).fillna(0)
    cathode_cross_trunk_weight[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_E_reco)\
          .astype(float).fillna(0)

    # UNCERTAINTY
    cathode_cross_trunk_weight_err = pd.Series(0, evtdf.index)
    cathode_cross_trunk_weight_err[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & ~cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & ~cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_W_split)*cathode_correction_err_W)\
          .astype(float).fillna(0)
    cathode_cross_trunk_weight_err[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & ~cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & ~cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_W_reco)*cathode_correction_err_W)\
          .astype(float).fillna(0)
    cathode_cross_trunk_weight_err[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & ~reco_trunk_crosses_cathode & cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_E_split)*cathode_correction_err_E)\
          .astype(float).fillna(0)
    cathode_cross_trunk_weight_err[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & cryoE] =\
      pd.cut(cathode_correction_trunk_var[true_trunk_crosses_cathode & reco_trunk_crosses_cathode & cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_E_reco)*cathode_correction_err_E)\
          .astype(float).fillna(0)

    # BRANCH
    true_branch_crosses_cathode = tpc_index(evtdf.branch.trk.truth.p.start.x) != tpc_index(evtdf.branch.trk.truth.p.end.x)
    true_same_cryo = (tpc_index(evtdf.branch.trk.truth.p.start.x)//2) == (tpc_index(evtdf.branch.trk.truth.p.end.x)//2)
    true_branch_crosses_cathode = true_branch_crosses_cathode & true_same_cryo
    reco_branch_crosses_cathode = tpc_index(evtdf.branch.trk.start.x) != tpc_index(evtdf.branch.trk.end.x)
    cryoE = evtdf.branch.trk.producer == 0
    
    # CORRECTION
    cathode_cross_branch_weight = pd.Series(1, evtdf.index)
    cathode_cross_branch_weight[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & ~cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & ~cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_W_split)\
          .astype(float).fillna(0)
    cathode_cross_branch_weight[true_branch_crosses_cathode & reco_branch_crosses_cathode & ~cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & reco_branch_crosses_cathode & ~cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_W_reco)\
          .astype(float).fillna(0)
    cathode_cross_branch_weight[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_E_split)\
          .astype(float).fillna(0)
    cathode_cross_branch_weight[true_branch_crosses_cathode & reco_branch_crosses_cathode & cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & reco_branch_crosses_cathode & cryoE], cathode_correction_bins, ordered=False, labels=cathode_correction_E_reco)\
          .astype(float).fillna(0)

    # UNCERTAINTY
    cathode_cross_branch_weight_err = pd.Series(0, evtdf.index)
    cathode_cross_branch_weight_err[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & ~cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & ~cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_W_split)*cathode_correction_err_W)\
          .astype(float).fillna(0)
    cathode_cross_branch_weight_err[true_branch_crosses_cathode & reco_branch_crosses_cathode & ~cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & reco_branch_crosses_cathode & ~cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_W_reco)*cathode_correction_err_W)\
          .astype(float).fillna(0)
    cathode_cross_branch_weight_err[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & ~reco_branch_crosses_cathode & cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_E_split)*cathode_correction_err_E)\
          .astype(float).fillna(0)
    cathode_cross_branch_weight_err[true_branch_crosses_cathode & reco_branch_crosses_cathode & cryoE] =\
       pd.cut(cathode_correction_branch_var[true_branch_crosses_cathode & reco_branch_crosses_cathode & cryoE], cathode_correction_bins, 
               ordered=False, labels=np.abs(1-cathode_correction_E_reco)*cathode_correction_err_E)\
          .astype(float).fillna(0)

    # Combine
    cathode_cross_weight = cathode_cross_trunk_weight * cathode_cross_branch_weight
    cathode_cross_weight_err = np.sqrt(cathode_cross_trunk_weight_err**2 +  cathode_cross_branch_weight_err**2)

    # Z=0 INDUCTION GAP
    gap_correction_trunk_var = evtdf.trunk.trk.dir.z
    gap_correction_branch_var = evtdf.branch.trk.dir.z
    
    # TRUNK
    true_trunk_crosses_gap = (evtdf.trunk.trk.truth.p.start.z >= 0) != (evtdf.trunk.trk.truth.p.end.z >= 0)
    reco_trunk_crosses_gap = (evtdf.trunk.trk.start.z >= 0) != (evtdf.trunk.trk.end.z >= 0)
    
    # CORRECTION
    gap_cross_trunk_weight = pd.Series(1, evtdf.index)
    gap_cross_trunk_weight[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & ~cryoE], gap_correction_bins, ordered=False, labels=gap_correction_W_split)\
          .astype(float).fillna(0)
    gap_cross_trunk_weight[true_trunk_crosses_gap & reco_trunk_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & reco_trunk_crosses_gap & ~cryoE], gap_correction_bins, ordered=False, labels=gap_correction_W_reco)\
          .astype(float).fillna(0)
    gap_cross_trunk_weight[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & cryoE], gap_correction_bins, ordered=False, labels=gap_correction_E_split)\
          .astype(float).fillna(0)
    gap_cross_trunk_weight[true_trunk_crosses_gap & reco_trunk_crosses_gap & cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & reco_trunk_crosses_gap & cryoE], gap_correction_bins, ordered=False, labels=gap_correction_E_reco)\
          .astype(float).fillna(0)

    # UNCERTAINTY
    gap_cross_trunk_weight_err = pd.Series(0, evtdf.index)
    gap_cross_trunk_weight_err[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & ~cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_W_split)*gap_correction_err)\
          .astype(float).fillna(0)
    gap_cross_trunk_weight_err[true_trunk_crosses_gap & reco_trunk_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & reco_trunk_crosses_gap & ~cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_W_reco)*gap_correction_err)\
          .astype(float).fillna(0)
    gap_cross_trunk_weight_err[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & ~reco_trunk_crosses_gap & cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_E_split)*gap_correction_err)\
          .astype(float).fillna(0)
    gap_cross_trunk_weight_err[true_trunk_crosses_gap & reco_trunk_crosses_gap & cryoE] =\
      pd.cut(gap_correction_trunk_var[true_trunk_crosses_gap & reco_trunk_crosses_gap & cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_E_reco)*gap_correction_err)\
          .astype(float).fillna(0)

    # BRANCH
    true_branch_crosses_gap = (evtdf.branch.trk.truth.p.start.z >= 0) != (evtdf.branch.trk.truth.p.end.z >= 0)
    reco_branch_crosses_gap = (evtdf.branch.trk.start.z >= 0) != (evtdf.branch.trk.end.z >= 0)
    
    # CORRECTION
    gap_cross_branch_weight = pd.Series(1, evtdf.index)
    gap_cross_branch_weight[true_branch_crosses_gap & ~reco_branch_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & ~reco_branch_crosses_gap & ~cryoE], gap_correction_bins, ordered=False, labels=gap_correction_W_split)\
          .astype(float).fillna(0)
    gap_cross_branch_weight[true_branch_crosses_gap & reco_branch_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & reco_branch_crosses_gap & ~cryoE], gap_correction_bins, ordered=False, labels=gap_correction_W_reco)\
          .astype(float).fillna(0)
    gap_cross_branch_weight[true_branch_crosses_gap & ~reco_branch_crosses_gap & cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & ~reco_branch_crosses_gap & cryoE], gap_correction_bins, ordered=False, labels=gap_correction_E_split)\
          .astype(float).fillna(0)
    gap_cross_branch_weight[true_branch_crosses_gap & reco_branch_crosses_gap & cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & reco_branch_crosses_gap & cryoE], gap_correction_bins, ordered=False, labels=gap_correction_E_reco)\
          .astype(float).fillna(0)

    # UNCERTAINTY
    gap_cross_branch_weight_err = pd.Series(0, evtdf.index)
    gap_cross_branch_weight_err[true_branch_crosses_gap & ~reco_branch_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & ~reco_branch_crosses_gap & ~cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_W_split)*gap_correction_err)\
          .astype(float).fillna(0)
    gap_cross_branch_weight_err[true_branch_crosses_gap & reco_branch_crosses_gap & ~cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & reco_branch_crosses_gap & ~cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_W_reco)*gap_correction_err)\
          .astype(float).fillna(0)
    gap_cross_branch_weight_err[true_branch_crosses_gap & ~reco_branch_crosses_gap & cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & ~reco_branch_crosses_gap & cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_E_split)*gap_correction_err)\
          .astype(float).fillna(0)
    gap_cross_branch_weight_err[true_branch_crosses_gap & reco_branch_crosses_gap & cryoE] =\
      pd.cut(gap_correction_branch_var[true_branch_crosses_gap & reco_branch_crosses_gap & cryoE], gap_correction_bins, 
               ordered=False, labels=np.abs(1-gap_correction_E_reco)*gap_correction_err)\
          .astype(float).fillna(0)

    # Combine
    gap_cross_weight = gap_cross_trunk_weight * gap_cross_branch_weight
    gap_cross_weight_err = np.sqrt(gap_cross_trunk_weight_err**2 +  gap_cross_branch_weight_err**2)

    return cathode_cross_weight, cathode_cross_weight_err, gap_cross_weight, gap_cross_weight_err

