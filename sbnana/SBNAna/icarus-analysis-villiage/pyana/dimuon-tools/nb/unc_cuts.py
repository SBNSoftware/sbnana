# Where Cuts are defined for the uncontained muons analysis.

# ----------------------------------------------------------
import pandas as pd 
import numpy as np 
from unc_funcs import *

# ----------------------------------------------------------
# Thresholds to use for Cuts

NuMI_angle_thresh = 3
open_angle_thresh = 15
max_shw_len_thresh = 0.

trk_len_thresh = 140. # cm

stub_dedx_l0_5cm_thresh = 50
stub_dedx_l1cm_thresh = 35
stub_dedx_l2cm_thresh = 15
stub_dedx_l3cm_thresh = 10

chi2_mu_forMuons = 13
chi2_p_forMuons = 90 

printstuff = False
if printstuff == True:
    print('\n')
    print('The default cut thresholds are:')
    print('NuMI_angle_thresh: ', NuMI_angle_thresh,' degrees')
    print('open_angle_thresh: ', open_angle_thresh,' degrees')
    print('max_shw_len_thresh: ', max_shw_len_thresh,' cm')
    print('trk_len_thresh: ', trk_len_thresh,' cm')
    print('stub_dedx_l0_5cm_thresh: ', stub_dedx_l0_5cm_thresh,' ')
    print('stub_dedx_l1cm_thresh: ', stub_dedx_l1cm_thresh,' ')
    print('stub_dedx_l2cm_thresh: ', stub_dedx_l2cm_thresh,' ')
    print('stub_dedx_l3cm_thresh: ', stub_dedx_l3cm_thresh,' ')
    print('chi2_mu_forMuons: ', chi2_mu_forMuons,' ')
    print('chi2_p_forMuons: ', chi2_p_forMuons,' ')
    print('\n')

# ---------------------------
# Cut Definitions
# Note: all these 'cuts' are actually just masks. The 'cut_all' function is what actually applies the cut.

# PID:

def is_muon(trk):
    return (trk.chi2pid.I2.chi2_muon < chi2_mu_forMuons) & (trk.chi2pid.I2.chi2_proton > chi2_p_forMuons)
def both_muon_tracks_mask(df, flip=False):
    mask = is_muon(df.trunk.trk) & is_muon(df.branch.trk) # this picks out the events where both are muons
    if flip==True:
        return ~mask, 'One or both candidate tracks are not muons'
    else:
        return mask, 'Both trunk, branch are muons'                                                        
## So that I can optimize the choice of chi2 values, rework the above cut into two separate ones:
##    (Using both of these cuts is the same as using both_muon_tracks_mask)
def ok_chi2mu(df, thresh = chi2_mu_forMuons):
    return( (df.trunk.trk.chi2pid.I2.chi2_muon < thresh) & (df.branch.trk.chi2pid.I2.chi2_muon < thresh) ), 'both chi2mu < %a' % chi2_mu_forMuons
def ok_chi2p(df, thresh = chi2_p_forMuons):
    return( (df.trunk.trk.chi2pid.I2.chi2_proton > thresh) & (df.branch.trk.chi2pid.I2.chi2_proton > thresh) ), 'both chi2p > %a' % chi2_p_forMuons

# KINEMATICS:

def numi_angle_mask(df, thresh = NuMI_angle_thresh, flip=False): # provide thresh in deg, but df has radians
    mask = df.Snumi_angle_wgtByLen < thresh*math.pi/180.
    if flip:
        return ~mask, 'S_NuMI_angle >= '+str(thresh)+' deg'#'\u00B0'
    else: 
        return mask, 'S_NuMI_angle < '+str(thresh)+' deg'#'\u00B0'

def open_angle_mask(df, thresh = open_angle_thresh, flip=False): # provide thresh in deg, but df has radians
    mask = np.arccos(dotdf(df.trunk.trk.dir, df.branch.trk.dir)) < thresh*math.pi/180.
    if flip:
        return ~mask, 'opening angle >= '+str(thresh)+' deg'
    else: 
        return mask, 'opening angle < '+str(thresh)+' deg'

# STUB ID:

def stub_mask(df, flip=False):
    hasstub = ((df.stub.l0_5cm.dqdx > stub_dqdx_l0_5cm_thresh) | 
               (df.stub.l1cm.dqdx > stub_dqdx_l1cm_thresh) | 
               (df.stub.l2cm.dqdx > stub_dqdx_l2cm_thresh) | 
               (df.stub.l3cm.dqdx > stub_dqdx_l3cm_thresh) 
              )
    if flip:
        return(hasstub), 'DOES have a stub (use for sideband)'
    else:
        return ~(hasstub), 'does not have stub'
# Recast into separate cuts for optimiztion:
def not_stub_05(df, thresh = stub_dedx_l0_5cm_thresh):
    hasstub = (df.stub.l0_5cm.dedx > thresh)
    return ~(hasstub), 'dE/dx <= %a MeV/cm up to 0_5 cm' % stub_dedx_l0_5cm_thresh
def not_stub_1(df, thresh = stub_dedx_l1cm_thresh):
    hasstub = (df.stub.l1cm.dedx > thresh)
    return ~(hasstub), 'dE/dx <= %a MeV/cm 0_5-1 cm' % stub_dedx_l1cm_thresh
def not_stub_2(df, thresh = stub_dedx_l2cm_thresh):
    hasstub = (df.stub.l2cm.dedx > thresh)
    return ~(hasstub), 'dE/dx <= %a MeV/cm 1-2 cm' % stub_dedx_l2cm_thresh
def not_stub_3(df, thresh = stub_dedx_l3cm_thresh):
    hasstub = (df.stub.l3cm.dedx > thresh)
    return ~(hasstub), 'dE/dx <= %a MeV/cm 2-3 cm' % stub_dedx_l3cm_thresh

# OBJECT CUTS:

def max_shw_len_mask(df, thresh = max_shw_len_thresh):
    #var = df.max_shw_len.copy()
    #var[np.isnan(var) | (var < 0)] = -10 
    has_long_shw = df.max_shw_len >= thresh
    return ~has_long_shw, 'max shower len < %a cm' % thresh

def longTrk_len_mask(df, thresh = trk_len_thresh):
    return var > thresh, 'longer of tracks > %a cm' % thresh
    
def shortTrk_len_mask(df, thresh = trk_len_thresh):
    return df.shorter_track_length > thresh, 'shorter of tracks > %a cm' % thresh

#def other_trk_len_mask(df, thresh = other_trk_len_thresh):
#    has_long_other_trk = df.max_othr_trk_len >= thresh
#    return ~has_long_other_trk, 'max other trk len < %a cm' % thresh

#def third_trk_dist_mask(df, thresh = third_trk_dist_thresh):
#    dontwant = df.third_trk_dist >= third_trk_dist_thresh
#    return ~dontwant, 'third_trk_dist < %a cm' % thresh

# -----------------------

# FUNCTION TO APPLY THE CUTS DEFINED ABOVE:

def apply_cuts(df, cuts, thresholds=None, sneaky_BG_cats=False, detailed_bsm=False, detailed_nu='none', hps_final_state=False, 
               flip_last_cut=False):
    
    #new_df = df.copy()
    if sneaky_BG_cats: categories = make_sneakyBG_categories(df)
    else: categories = make_categories(df, detailed_bsm=detailed_bsm, detailed_nu=detailed_nu, hps_final_state=hps_final_state)
    
    # initialze data frames
    cut_results_df = pd.DataFrame(
        np.zeros((0,len(categories)), dtype=int), # start w/ zero rows (cuts), fill later
        columns = [c.name for c in categories]
    )
    cut_results_df_mc = cut_results_df.copy() #deep=True
    cut_results_df_pot = cut_results_df.copy() #deep=True
    cut_results_df_percent = cut_results_df.copy() #deep=True
    
    # fill in first row of data frame for "no cuts"
    row_mc = []
    row_pot = []
    for c in categories:
        #print(sum(df[c].scale))
        row_mc.append(df[c].shape[0])
        row_pot.append(round(100*sum(df[c].scale))/100.)
    cut_results_df_mc.loc["preselection"] = row_mc 
    first_row_mc = row_mc
    cut_results_df_pot.loc["preselection"] = row_pot
    cut_results_df_percent.loc["preselection"] = [1.] * len(row_pot)
    
    # Loop through cuts to make rows for data frame and to make master_mask
    master_mask = None
    for i in range(len(cuts)):
        if thresholds is None:
            func_output = cuts[i](*[df])
        else: 
            func_output = cuts[i](*[df], thresh=thresholds[i])
        if flip_last_cut: #overwrite with the flip
            if i == len(cuts)-1:
                func_output = cuts[i](*[df], flip=True)
        if i==0: 
            master_mask = func_output[0]
        else:
            master_mask = master_mask & func_output[0]
        new_df = df[master_mask]
        if sneaky_BG_cats: new_categories = make_sneakyBG_categories(new_df)
        else: new_categories = make_categories(new_df, detailed_bsm=detailed_bsm, detailed_nu=detailed_nu, hps_final_state=hps_final_state)
        
        row_mc = []
        row_pot = []
        for c in new_categories:
            row_mc.append(new_df[c].shape[0])
            row_pot.append(round(100*sum(new_df[c].scale))/100.)
        cut_results_df_mc.loc[func_output[1]] = row_mc
        cut_results_df_pot.loc[func_output[1]] = row_pot
        cut_results_df_percent.loc[func_output[1]] = np.array(row_mc)/np.array(first_row_mc)
    
    return cut_results_df_mc, cut_results_df_pot, cut_results_df_percent, master_mask 

# -----------------------

# DICTIONARY OF DIFFERENT EVENT SELECTIONS (CUTS + THRESHOLDS):
# Used when running run_evtSel_script.py
# Instructions: It is important to NEVER REPLACE OR DELETE selections in this dictionary that you've already run.

evtSel_dict = dict([
    ( "A",# The 2024 benchmarked evtsel.
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu, ok_chi2p,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            0.0,
            13, 90,
            3, 15,
            140.0
        ]
     )
    ),

    ( "B", # Bigger opening angle
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu, ok_chi2p,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            0.0,
            13, 90,
            3, 25,
            140.0
        ]
     )
    ),
    
    ( "C", # Bigger opening angle; Longer trk len
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu, ok_chi2p,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            0.0,
            13, 90,
            3, 25,
            170.0
        ]
     )
    ), 
    
    ( "D",# Loosen theta_NuMI
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu, ok_chi2p,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            0.0,
            13, 90,
            5, 15,
            140.0
        ]
     )
    ),
    
    ( "E",# Loosen theta_NuMI, Longer trk len
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu, ok_chi2p,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            0.0,
            13, 90,
            5, 15,
            170.0
        ]
     )
    ),
    
    ( "F",# Remove Shw Len Cut
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            ok_chi2mu, ok_chi2p,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            13, 90,
            3, 15,
            140.0
        ]
     )
    ),
    
    ( "G",# Tighten Shw Len Cut
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu, ok_chi2p,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            10.0,
            13, 90,
            3, 15,
            140.0
        ]
     )
    ),
    
    ( "H",# Remove Chis_proton Cut
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            0.0,
            13,
            3, 15,
            140.0
        ]
     )
    ), 
    
    ( "I",# Combine Selections G and H: loosen shower length cut and drop chis_proton cut.
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu,
            numi_angle_mask, open_angle_mask,
            shortTrk_len_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            10.0,
            13,
            3, 15,
            140.0
        ]
     )
    ),
    
    ###
    # May 19, 2025
    # Everything above were used in studies I performed Jan/Feb 2025. Don't touch them!
    # Below here I am playing around with potential near sidebands.
    ###
    
    # check mc sideband in theta_NuMI:
    ( "J",
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu,
            open_angle_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            10.0,
            13,
            15
        ]
     )
    ),
    
    # check mc sideband in theta_mumu:
    ( "K",
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu,
            numi_angle_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            10.0,
            13,
            3
        ]
     )
    ),
    
    # check mc sideband in theta_mumu:
    ( "L",
     ( # tuple of cut list and thresholds
        [ # cut list
            not_stub_05, not_stub_1, not_stub_2, not_stub_3,
            max_shw_len_mask,
            ok_chi2mu,
            numi_angle_mask, open_angle_mask
        ],
        [ # thresholds
            50, 35, 15, 10,
            10.0,
            13,
            3, 15
        ]
     )
    )

])
                           