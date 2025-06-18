# functions copied from unc_MC Event Selection.ipynb on September 28, 2023

import math
from math import dist
import numpy as np
import pandas as pd 
from util import *
from pyanalib import panda_helpers
from numpy import sqrt, cos, sin
import var

GOAL_POT = 2.41e20 # combined R1 and R2.
POTSTR = "2.41e20 POT"

nmu = var.DF.iscc & ((var.DF.pdg == 14) | (var.DF.pdg == -14))
npi = var.DF.npi
ns = var.DF.nsm + var.DF.nsp
is_coh_like = (nmu + npi + ns >= 2) & (var.DF.max_proton_ke < 0.05)

def is_coh_like_JD(df):
    nmu = df.slc.truth.nmu # df.slc.truth.iscc & ((df.slc.truth.pdg == 14) | (df.slc.truth.pdg == -14))
    npi = df.slc.truth.npi
    ns = df.slc.truth.nsm + df.slc.truth.nsp # Sigma+ or Sigma-
    #print(type(df.slc.truth.max_proton_ke < 0.05))
    return (nmu + npi + ns >= 2) & (df.slc.truth.max_proton_ke < 0.05)

def is_NOT_coh_like_JD(df):
    return ~is_coh_like_JD(df)

def satisfies_new_FV(df): # df should be evtdf
    return ( 
        SlcInFV(df.slc.vertex) & 
        ( (TrkInFV(df.trunk.trk.end)) | (df.trunk.trk.len >= 100)) & 
        ( (TrkInFV(df.branch.trk.end)) | (df.branch.trk.len >= 100)) 
    ) # each track is either contained, or long enough to be uncontained but still kept.
# July 22, 2024. This is needed because the FV definition has changed since the preselection was run on the cafs to make the dataframes. It is now more stringent, so there are some slices (tracks) in the dataframes that should not have been selected given the updated FV definition.

def decay_in_icarus(in_dist, out_dist, mean_dist): # This takes care of the integral in the decay weight. # same as what Gray defined as "decay weight." I've changed the name bc technically the decay weight also includes a factor of the branching fraction (which is just 1 for hps, where it always decays into muons where we are looking).
    return np.exp(-in_dist / mean_dist) - \
          np.exp(-out_dist / mean_dist)

def flux_weight(mixing): 
    return mixing * mixing

def reweight_mixing(newmixing, start, enter, exit, mean_dist, oldmixing=1e-5):
    # Use this one for HPS. 
    #Technically, only scaling 1/fa in the ALPs model also follows this form. But in general, cl can scale too. So use next function when reweighting ALPs.
    dist_in = dist(enter, start)
    dist_out = dist(exit, start)
    
    old_dcy = decay_in_icarus(dist_in, dist_out, mean_dist)
    new_dcy = decay_in_icarus(dist_in, dist_out, mean_dist / (newmixing**2 / oldmixing**2))
    
    old_flux = flux_weight(oldmixing)
    new_flux = flux_weight(newmixing)
    
    return (new_flux / old_flux) * (new_dcy / old_dcy)

def reweight_alps(old_fa, new_fa, oldcl, newcl, start, enter, exit, mean_dist, oldf, print_stuff=False): # f is branching fraction.
# Use this function to rescale generated alps to different alp benchmarks.
    dist_in = dist(enter, start)
    dist_out = dist(exit, start)
    old_decay_weight = decay_in_icarus(dist_in, dist_out, mean_dist)*oldf
    new_mean_dist = mean_dist*(new_fa**2/old_fa**2)/( (1-oldf) + (oldf*newcl**2/oldcl**2) )
    new_f = (newcl**2/oldcl**2)*oldf/( (newcl**2/oldcl**2)*oldf + 1-oldf )
    new_decay_weight = decay_in_icarus(dist_in, dist_out, new_mean_dist)*new_f
    
    flux_weight_rescaling = old_fa**2/new_fa**2  

    if print_stuff:
        print('cl, fa: ', newcl, ', ', new_fa)
        print('decay_weight: ', float(new_decay_weight))
        print('mean_dist: ', float(new_mean_dist))
        print('decay in icarus: ', float(decay_in_icarus(dist_in, dist_out, new_mean_dist)))
        print('f: ', float(new_f))
    
    return flux_weight_rescaling*new_decay_weight/old_decay_weight

# FUNCTION TO REWEIGHT HPSs TO ALPs.

def reweight_hps_to_alps(mass, hps_mix, hps_mean_dist, hps_f, 
                         alp_fa, alp_cl, alp_f, # ("new" alp_f, cause there is no "old" alp_f when starting with hps)
                         start, enter, exit, kaon_parent_pdg,
                         alp_c1=1, alp_c2=1, alp_c3=1,
                         k_to_alp_strong = True
                        ): # f's are branching fractions to muons relative to total decay width for hps or alp.
    dist_in = dist(enter, start)
    dist_out = dist(exit, start)
    
    # Starting up tomorrow (Thursday 10/31/24): fill in alp_flux_weight_piece - strong or weak? Answer: strong!
    
    if kaon_parent_pdg == 130: # K0L
        hps_flux_weight_piece = HPS_KaonLongBranchingRatio(mass, hps_mix) 
        if k_to_alp_strong:
            alp_flux_weight_piece =  ALP_KLongStrongBranchingRatio(mass, alp_fa, alp_c3) # strong
        else:
            alp_flux_weight_piece = ALP_KLongWeakBranchingRatio(mass, alp_fa, alp_c2) # weak
    elif abs(kaon_parent_pdg) == 321: # K+/-
        hps_flux_weight_piece = HPS_KaonPlusBranchingRatio(mass, hps_mix) 
        if k_to_alp_strong:
            alp_flux_weight_piece = ALP_KPlusStrongBranchingRatio(mass, alp_fa, alp_c3) # strong
        else:
            alp_flux_weight_piece = ALP_KPlusWeakBranchingRatio(mass, alp_fa, alp_c2) # weak
    else: 
        print("We've got a problem here!")
    hps_decay_weight_piece = decay_in_icarus(dist_in, dist_out, hps_mean_dist)*hps_f
    alp_mean_dist = hps_mean_dist*HPS_total_decay_Width(mass, hps_mix)/ALP_total_decay_Width(mass, alp_fa, alp_cl, alp_c1, alp_c2, alp_c3)
    alp_decay_weight_piece = decay_in_icarus(dist_in, dist_out, alp_mean_dist)*alp_f
    
    return (alp_flux_weight_piece/hps_flux_weight_piece)*(alp_decay_weight_piece/hps_decay_weight_piece)
   
    
def jamie_sample_concat(dflist): # dflist should be a list of multi-index dataframes
    master_df = pd.concat(dflist, keys=np.arange(len(dflist)), names=['sample'])
    return master_df

# The following function takes an event dataframe, and uses truth information to determine what type of
#    neutrino interaction it is. It then saves that as 'nu_mode' and adds the label to the dataframe, returning the df.
# WARNING: only apply this to a dataframe with only neutrinos, otherwise you might get confused.
def add_nu_int_type(df):
    new_df = df.copy()
    nu_mode = []
    for ind in new_df.index:
        row = df.loc[ind]
        iscc = bool(int(row.slc.truth.iscc))
        genie_mode = int(row.slc.truth.genie_mode)
        nu_pdg = int(row.slc.truth.pdg)
        npi = int(row.slc.truth.npi)
        npi0 = int(row.slc.truth.npi0)
        max_p_ke = int(row.slc.truth.max_proton_ke)
        if (np.abs(nu_pdg) == 14) & (iscc): # nu mu CC
            if genie_mode == 3: # nu mu CC COH 
                nu_mode.append('$\\nu_\\mu$ CC COH')
            elif npi >=1 & npi0 == 0:
                if  max_p_ke < 0.02:
                    nu_mode.append("$\\nu_\\mu$ CC n$\\pi$0p")
                else:
                    nu_mode.append("$\\nu_\\mu$ CC n$\\pi$np")
            else: 
                nu_mode.append("$\\nu_\\mu$ CC Other") 
        else: 
            nu_mode.append("not $\\nu_\\mu$ CC") 
        #except: nu_mode.append('-')
    new_df['nu_mode'] = nu_mode
    return new_df
    
def add_hdr_info(df, hdrs):
    new_df = df.copy()
    proc = []
    clust = []
    for idx in new_df.index:
        p = hdrs[idx[0]].loc[(idx[1],idx[2])].proc
        c = hdrs[idx[0]].loc[(idx[1],idx[2])].cluster
        proc.append(p)
        clust.append(c)
    new_df['proc'] = proc
    new_df['cluster'] = clust
    return new_df

def simpler_add_hdr_info(df, hdr): # for a single df and its header (not conatenated df and list of hdrs like the other one is)
    new_df = df.copy()
    proc = []
    clust = []
    daqrun = []
    evt_from_hdr = []
    for idx in new_df.index:
        proc.append(hdr.loc[(idx[0],idx[1])].proc)
        clust.append(hdr.loc[(idx[0],idx[1])].cluster)
        daqrun.append(hdr.loc[(idx[0],idx[1])].run)
        evt_from_hdr.append(hdr.loc[(idx[0],idx[1])].evt)
    new_df['proc'] = proc
    new_df['cluster'] = clust
    new_df['daqrun'] = daqrun
    new_df['evt_from_hdr'] = evt_from_hdr
    return new_df

# COLORS
# Good tools for picking colors: 
# https://imagecolorpicker.com
# https://htmlcolorcodes.co# Use this to find html hex string color codes: https://htmlcolorcodes.com 


blues = ["#B0E0E6", "#87CEEB", "#6495ED", "#1E90FF","#4682B4", "#00008B", "#0000FF", "#4169E1", "#1E90FF", "#4682B4"]
greens = ["#CFFCDA", "#8DF9A7", "#72DF8C", "#48B161", "#2C8942", "#145C25", "#054915", "#D1F725", "#00FF00", "#32CD32", "#3CB371", "#008000"]
purples = ["#D8BFD8", "#9370DB", "#663399"]
oranges = ["#FECC88", "#FCC171", "#ECA646", "#DE9025", "#CD7A08", "#A15F03"]#, "#", "#"]
#greens = ["#8FBC8F", "#2E8B57", "#556B2F", "#006400"]
#oranges = ["#FFDAB9", "#A0522D", "#F4A460"]

GraysColors = ["#EAC387", "#97B9BB", "#C6C2D9", "#F0CBD2", "#983E49", "#A57647", "lightgray"]
southwest = ["#caa834", "#1e7195", "#ff8942", "#d2998f", "#4aede4", "#2b6751", "#946671"]
soladero_lime = ["#caa834", "#e03b14", "#1e7195", "#d2998f", "#3d6f68", "#7ddc69", "#946671"]
soladero = ["#caa834", "#e03b14", "#1e7195", "#d2998f", "#3d6f68", "#adceb4", "#946671"]

# BROAD CATEGORIES ( BSM model + benchmark | nu | cosmic )
# tmatch = truth match (based on what contributed most E to slice)

# copied from event selection nb 8/29/24, and also used for sensitivities notebook.
def make_categories(df, detailed_bsm=False, detailed_nu='none', hps_final_state=False): 
# Use this function to categorize the MC samples (signal and bg) that are used in designing the event selection.
# hps_final_state is used to distinguish for mumu and pipi events. 
    
    # COSMICS
    
    is_cosmic = df.slc.tmatch.idx < 0 # from any file, so I'm assuming both samples include cosmics.
    is_cosmic.name = "Cosmic"
    is_cosmic.color = soladero_lime[6] #"C3"
    
    # SIGNAL
    
    if detailed_bsm:
        
        is_higgs = ( (df.slc.tmatch.idx >= 0) & df.higgs & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0) )# Only consider muon channel, exclude any Higgs that decayed to pions.
        #is_higgs.name = "Scalar"
        higgs_benchmarks = []
        cat_hdrs = []
        #higgs_benchmark_names = df[is_higgs]["sample"].unique()
        higgs_benchmark_names = df[is_higgs]["sample"].unique() # Made this change June 3 to avoid bug when applying cuts, where you eliminate all of one category so your results dataframe shape gets messed up. Hopefully this doesn't cause other problems.
        for l in range(len(higgs_benchmark_names)):
            if hps_final_state:
                bm_mu = is_higgs & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0) & (df["sample"] == higgs_benchmark_names[l])
                bm_mu.name = higgs_benchmark_names[l] + ' $\\mu\\mu$'
                bm_mu.color = blues[l]
                higgs_benchmarks.append(bm_mu)
                if float(df[is_higgs & (df["sample"] == higgs_benchmark_names[l])].iloc[[0]]["bsm_mass"] > (2*135.)):
                    bm_pipm = is_higgs & (df.slc.truth.nmu == 0) & (df.slc.truth.npi0 == 0) & (df["sample"] == higgs_benchmark_names[l])
                    bm_pipm.name = higgs_benchmark_names[l] + ' $\\pi\\pi +-$'
                    bm_pipm.color = blues[l]
                    higgs_benchmarks.append(bm_pipm)
                    bm_pi0 = is_higgs & (df.slc.truth.nmu == 0) & (df.slc.truth.npi == 0) & (df["sample"] == higgs_benchmark_names[l])
                    bm_pi0.name = higgs_benchmark_names[l] + ' $\\pi\\pi 0$'
                    bm_pi0.color = blues[l]
                    higgs_benchmarks.append(bm_pi0)
            else:
                bm = is_higgs & (df["sample"] == higgs_benchmark_names[l])
                bm.name = higgs_benchmark_names[l]
                bm.color = blues[l]
                higgs_benchmarks.append( bm )    
        
        is_alp = (df.slc.tmatch.idx >= 0) & df.alp & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0)
        #is_alp.name = "ALP"
        alp_benchmarks = []
        alp_benchmark_names = df[is_alp]["sample"].unique()
        for l in range(len(alp_benchmark_names)):
            bm = is_alp & (df["sample"] == alp_benchmark_names[l])
            bm.name = alp_benchmark_names[l]# + " (no sup.)"
            bm.color = greens[l]
            alp_benchmarks.append( bm )
        
        bsm_cats = higgs_benchmarks + alp_benchmarks
            
    else:
        is_bsm = (df.slc.tmatch.idx >= 0) & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0) & ( df.higgs | df.alp ) # ( df.higgs | df.alp_withsup | df.alp_nosup )
        is_bsm.name = "BSM"
        is_bsm.color = "#191970" #, "#00FF7F" #"#000000"
        bsm_cats = [is_bsm]

      
    # NEUTRINOS
    
    is_nu = (df.slc.tmatch.idx >= 0) & df.nu

    if detailed_nu == 'int_type':
        
        nu_NC = is_nu & ((df.slc.truth.iscc == 0) & (df.slc.tmatch.idx >= 0))
        nu_NC.name = "$\\nu$ NC"
        nu_NC.color = soladero_lime[0] #"#caa834" #"#EAC387"
        
        numu_CC_QE_MEC = is_nu & ((df.slc.truth.pdg == 14) | (df.slc.truth.pdg == -14)) & (df.slc.truth.iscc==1) & ((df.slc.truth.genie_mode == 0) | (df.slc.truth.genie_mode == 10)) # CC QE+MEC
        numu_CC_QE_MEC.name = "$\\nu_\\mu$ CC QE+MEC"
        numu_CC_QE_MEC.color = soladero_lime[1] 
        
        numu_CC_RES = is_nu & ((df.slc.truth.pdg == 14) | (df.slc.truth.pdg == -14)) & (df.slc.truth.iscc==1) & (df.slc.truth.genie_mode == 1) # CC RES
        numu_CC_RES.name = "$\\nu_\\mu$ CC RES"
        numu_CC_RES.color = soladero_lime[2]
        
        numu_CC_DIS = is_nu & ((df.slc.truth.pdg == 14) | (df.slc.truth.pdg == -14)) & (df.slc.truth.iscc==1) & (df.slc.truth.genie_mode == 2) # CC DIS
        numu_CC_DIS.name = "$\\nu_\\mu$ CC DIS"
        numu_CC_DIS.color = soladero_lime[3] 
        
        numu_CC_COH = is_nu & ((df.slc.truth.pdg == 14) | (df.slc.truth.pdg == -14)) & (df.slc.truth.iscc==1) & (df.slc.truth.genie_mode == 3) # CC COH
        numu_CC_COH.name = "$\\nu_\\mu$ CC COH"
        numu_CC_COH.color = soladero_lime[4]
        
        nu_other = is_nu & ~nu_NC & ~numu_CC_QE_MEC & ~numu_CC_RES & ~numu_CC_DIS & ~numu_CC_COH
        nu_other.name = "$\\nu$ Other"
        nu_other.color = soladero_lime[5]
        
        nu_cats = [nu_NC] + [numu_CC_QE_MEC] + [numu_CC_RES] + [numu_CC_DIS] + [numu_CC_COH] + [nu_other] 
        
    elif detailed_nu == 'final_state':
        
        numu_cc_coh = is_nu & (df.slc.truth.genie_mode == 3) & (np.abs(df.slc.truth.pdg) == 14) & (df.slc.truth.iscc.astype('bool'))
        numu_cc_coh.name = '$\\nu_\\mu$ CC COH'
        numu_cc_coh.color = '#FADD28' #oranges[0]
        
        numu_cc_npizp = (is_nu & 
                         (df.slc.truth.genie_mode != 3) & (np.abs(df.slc.truth.pdg) == 14) & 
                         (df.slc.truth.iscc.astype('bool')) & (df.slc.truth.npi >= 1) & (df.slc.truth.npi0==0) &
                         (df.slc.truth.max_proton_ke < 0.02) )
        numu_cc_npizp.name = "$\\nu_\\mu$ CC n$\\pi$0p"
        numu_cc_npizp.color = oranges[1]
        
        numu_cc_npinp = ( is_nu & 
                         (df.slc.truth.genie_mode != 3) & (np.abs(df.slc.truth.pdg) == 14) & 
                         (df.slc.truth.iscc.astype('bool')) & (df.slc.truth.npi >= 1) & (df.slc.truth.npi0==0) & 
                         (df.slc.truth.max_proton_ke >= 0.02))
        numu_cc_npinp.name = '$\\nu_\\mu$ CC n$\\pi$np'
        numu_cc_npinp.color = oranges[2]
        
        numu_cc_other = ( is_nu & 
                         (df.slc.truth.genie_mode != 3) & 
                         (np.abs(df.slc.truth.pdg) == 14) & (df.slc.truth.iscc.astype('bool')) & 
                         ( (df.slc.truth.npi<1) | (df.slc.truth.npi0!=0) )
                        )
        numu_cc_other.name = '$\\nu_\\mu$ CC Other'
        numu_cc_other.color = oranges[3]
        
        not_numu_cc = is_nu & ( (np.abs(df.slc.truth.pdg) != 14) | ~df.slc.truth.iscc.astype('bool') )
        not_numu_cc.name = 'not $\\nu_\\mu$ CC'
        not_numu_cc.color = oranges[4]
        
        nu_cats = [numu_cc_coh] + [numu_cc_npizp] + [numu_cc_npinp] + [numu_cc_other] + [not_numu_cc]
        
    else:
        is_nu.name = "$\\nu$"
        is_nu.color = "C1"
        nu_cats = [is_nu]

    categories = bsm_cats + nu_cats + [is_cosmic]
    return categories

# copied from event selection nb 9/6/24
def make_vtxcategories(df):
    
    trunkIDcategories = [
        (df.trunk.trk.truth.p.pdg == 13) | (df.trunk.trk.truth.p.pdg == -13), # "$\\mu$"
        (df.trunk.trk.truth.p.pdg == 211) | (df.trunk.trk.truth.p.pdg == -211), # "$\\pi$"
        (df.trunk.trk.truth.p.pdg == 2212), # "$p$"
        (df.trunk.trk.truth.p.pdg == 22), # "$\\gamma$"
    ]
    trunkIDcategories.append(~trunkIDcategories[0] & ~trunkIDcategories[1] & ~trunkIDcategories[2] & ~trunkIDcategories[3])

    
    branchIDcategories = [
        (df.branch.trk.truth.p.pdg == 13) | (df.branch.trk.truth.p.pdg == -13), # "$\\mu$"
        (df.branch.trk.truth.p.pdg == 211) | (df.branch.trk.truth.p.pdg == -211), # "$\\pi$"
        (df.branch.trk.truth.p.pdg == 2212), # "$p$"
        (df.branch.trk.truth.p.pdg == 22),
    ]
    branchIDcategories.append(~branchIDcategories[0] & ~branchIDcategories[1] & ~branchIDcategories[2] & ~branchIDcategories[3])

    
    vtxIDcategories = [
        trunkIDcategories[0] & branchIDcategories[0], # mu mu
        (trunkIDcategories[0] & branchIDcategories[1]) | (trunkIDcategories[1] & branchIDcategories[0]), # mu pi+-
        (trunkIDcategories[0] & branchIDcategories[2]) | (trunkIDcategories[2] & branchIDcategories[0]), # mu p
        (trunkIDcategories[1] & branchIDcategories[2]) | (trunkIDcategories[2] & branchIDcategories[1]) # p pi+-
        # mu pi0
        # mu gamma
        # pi+- pi+-
    ]
    vtxIDcategories.append(~vtxIDcategories[0] & ~vtxIDcategories[1] & ~vtxIDcategories[2] & ~vtxIDcategories[3])
    
    vtxIDcategories[0].name = "$\\mu\\mu$"
    vtxIDcategories[1].name = "$\\mu\\pi$"
    vtxIDcategories[2].name = "$\\mu p$"
    vtxIDcategories[3].name = "$p\\pi$"
    vtxIDcategories[4].name = "other"
    
    vtxIDcategories[0].color = southwest[0]
    vtxIDcategories[1].color = southwest[1]
    vtxIDcategories[2].color = southwest[2]
    vtxIDcategories[3].color = '#CE3106' #southwest[3]
    vtxIDcategories[4].color = southwest[4]
    
    return vtxIDcategories

def make_sneakyBG_categories(df):
    # This funciton is based on Truth, and the categories are not orthogonal.
    # Categories here are designed to help me see what kinds of background different cuts target.
    
    trunkIDcategories = [
        ( abs(df.trunk.trk.truth.p.pdg) == 13 ), # "$\\mu$"
        ( abs(df.trunk.trk.truth.p.pdg) == 211), # "$\\pi+-$"
        (df.trunk.trk.truth.p.pdg == 2212), # "$p$"
        (df.trunk.trk.truth.p.pdg == 22), # "$\\gamma$" 
    ]
    trunkIDcategories.append(
        ~trunkIDcategories[0] & ~trunkIDcategories[1] & ~trunkIDcategories[2] & ~trunkIDcategories[3]
    ) # other

    
    branchIDcategories = [
        ( abs(df.branch.trk.truth.p.pdg) == 13 ), # "$\\mu$"
        ( abs(df.branch.trk.truth.p.pdg) == 211), # "$\\pi+-$"
        (df.branch.trk.truth.p.pdg == 2212), # "$p$"
        (df.branch.trk.truth.p.pdg == 22), # "$\\gamma$"
    ]
    branchIDcategories.append(
        ~branchIDcategories[0] & ~branchIDcategories[1] & ~branchIDcategories[2] & ~branchIDcategories[3]
    ) # other

    
    vtxIDcategories = [
        branchIDcategories[1] | trunkIDcategories[1], # has a pi+-
        branchIDcategories[2] | trunkIDcategories[2], # has a p
        branchIDcategories[3] | trunkIDcategories[3], # has a gamma
        branchIDcategories[4] | trunkIDcategories[4], # has something other than mu, pi+-, p, or gamma
        ( df.slc.truth.max_proton_ke > 0. ) & ( df.slc.truth.max_proton_ke < 0.5 ) # has a low E proton
    ]
    
    vtxIDcategories[0].name = "has a $\\pi+-$"
    vtxIDcategories[1].name = "has a $p$"
    vtxIDcategories[2].name = "has a $\\gamma$"
    vtxIDcategories[3].name = "has an exotic" #"has something not $\\mu$, $\\pi+-$, $p$, or $\\gamma$"
    vtxIDcategories[4].name = "has a low E $p$"
    
    #vtxIDcategories = [
    #    trunkIDcategories[0] & branchIDcategories[0], # mu mu
    #    (trunkIDcategories[0] & branchIDcategories[1]) | (trunkIDcategories[1] & branchIDcategories[0]), # mu pi+-
    #    (trunkIDcategories[0] & branchIDcategories[2]) | (trunkIDcategories[2] & branchIDcategories[0]), # mu p
    #    (trunkIDcategories[1] & branchIDcategories[2]) | (trunkIDcategories[2] & branchIDcategories[1]), # p pi+-
    #    # mu pi0
    #    (trunkIDcategories[3] & branchIDcategories[1]) | (trunkIDcategories[3] & branchIDcategories[0]) # mu gamma
    #    # pi+- pi+-
    #]
    #vtxIDcategories.append(~vtxIDcategories[0] & ~vtxIDcategories[1] & ~vtxIDcategories[2] & ~vtxIDcategories[3])
    
    #vtxIDcategories[0].name = "$\\mu\\mu$"
    #vtxIDcategories[1].name = "$\\mu\\pi$"
    #vtxIDcategories[2].name = "$\\mu p$"
    #vtxIDcategories[3].name = "$p\\pi$"
    #vtxIDcategories[4].name = "other"
    
    vtxIDcategories[0].color = southwest[0]
    vtxIDcategories[1].color = southwest[1]
    vtxIDcategories[2].color = southwest[2]
    vtxIDcategories[3].color = '#CE3106' #southwest[3]
    vtxIDcategories[4].color = southwest[4]
    
    return vtxIDcategories

def apply_cuts(df, cuts, detailed_hps=False, flip_last_cut=False):
    
    #new_df = df.copy()
    categories = make_categories(df, detailed_hps=detailed_hps)
    
    # initialze data frames
    cut_results_df = pd.DataFrame(
        np.zeros((0,len(categories)), dtype=int), # start w/ zero rows (cuts), fill later
        columns = [c.name for c in categories] # [c.name.split(',')[-1] for c in categories]
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
        row_pot.append(sum(df[c].scale))
    cut_results_df_mc.loc["preselection"] = row_mc 
    first_row_mc = row_mc
    cut_results_df_pot.loc["preselection"] = row_pot
    cut_results_df_percent.loc["preselection"] = [1.] * len(row_pot)
    
    # Loop through cuts to make rows for data frame and to make master_mask
    for i in range(len(cuts)):
        func_output = cuts[i](*[df])
        if flip_last_cut: #overwrite with the flip
            if i == len(cuts)-1:
                func_output = cuts[i](*[df], flip=True)
                
        if i==0: 
            master_mask = func_output[0]
        else:
            master_mask = master_mask & func_output[0]
        new_df = df[master_mask]
        new_categories = make_categories(new_df, detailed_hps=detailed_hps)
        
        row_mc = []
        row_pot = []
        for c in new_categories:
            row_mc.append(new_df[c].shape[0])
            row_pot.append(sum(new_df[c].scale))
        cut_results_df_mc.loc[func_output[1]] = row_mc   
        cut_results_df_pot.loc[func_output[1]] = row_pot
        cut_results_df_percent.loc[func_output[1]] = np.array(row_mc)/np.array(first_row_mc)
    
    return cut_results_df_mc, cut_results_df_pot, cut_results_df_percent, master_mask 

def distance_3d(point1, point2):
    return np.linalg.norm(point2 - point1)

def angle_between_vecs(a,b, deg=False): #a and b need to be numpy arrays of same length
    rad = np.arccos(a.dot(b)/(np.sqrt(a.dot(a))*np.sqrt(b.dot(b))))
    if deg:
        return rad*180./np.pi
    else:
        return rad

# NuMI Angle stuff:
# Geometry stuff is copied from here:
# https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/EventGenerator/MeVPrtl/config/numi_kaon_common.fcl#L7C15-L7C101
# Which is copied from here:
# https://github.com/SBNSoftware/icaruscode/blob/develop/fcl/gen/numi/GNuMIFlux.xml#L765
# The geometry is described in the TN:
# https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=27384
# Coordinate transformation from NuMI to ICARUS system (which I used to get NUMIDIR) is documented at:
# https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=22998&filename=main_v2.pdf&version=2 (this is referenced in TN)

beam2det = np.array([[0.921035925, 0.022715103, 0.388814672],
                     [0, 0.998297825, -0.058321970],
                     [-0.389477631, 0.053716629, 0.919468161]])
beamorigin = np.array([4.503730e2, 80.153901e2, 795.112945e2]) 
# This the vector from beam origin (NuMI target) to ICARUS detector (origin), in beam coordinates.

BEAMDIR = beam2det.dot(beamorigin) / np.linalg.norm(beam2det.dot(beamorigin)) 
# this is the vector from NuMI target to ICARUS detector, in ICARUS coordinates.
#print(BEAMDIR)

NUMIDIR = beam2det.dot(np.array([0,0,1.]))/np.linalg.norm(beam2det.dot(np.array([0,0,1.]))) 
# This is direction of NuMI beamline, in ICARUS coordinates.

absorber_origin = beamorigin - np.array([0, 0, 730e2]) 
# this is vector from NuMI Absorber to ICARUS detector, in beam coordinates
ABSORBDIR = beam2det.dot(absorber_origin) / np.linalg.norm(beam2det.dot(absorber_origin))
# this is vector from NuMI Absorber to ICARUS detector, in ICARUS coordinates.

def Sbeamangle(trunk_track, branch_track, beamdir, method): # This is what I call theta_NuMI # Use this for Reco!!
    if method == 'planar': 
        # method: (Note: this method works no matter what order you take the cross product - good!)
        #         (Note: this method introduces ambiguity between truly acute or truly obtuse angles wrt BEAMDIR, and related to that, the range in returned values is -90 to 90 and then I take abs so actually 0 to 90.)
        # 1) take cross product of track directions to get a vector normal to their plane.
        # 2) get angle between the normal vector and beamdir.
        # 3) angle between normal vector and beamdir is 90 deg. minus the above angle.
        trk_xprd = np.cross(trunk_track.dir, branch_track.dir) # cross product verified as working!
        trk_xprd_mag = np.linalg.norm(trk_xprd, axis=1)
        trk_xprd_norm = trk_xprd/trk_xprd_mag[:, None]
        angle1 = np.arccos(np.dot(trk_xprd_norm, BEAMDIR)) # these range from 0 to pi.
        Sbeamangle = abs(angle1 - math.pi/2) # these range from -pi/2 to pi/2, which is a different range than what the other func methods would give (0, pi)
        return Sbeamangle
    if method == 'p0_p1_truth': #input would be eg evtdf.slc.truth.p0
        scalar_mom = trunk_track.genp + branch_track.genp
    elif method == 'track_truth':
        scalar_mom = trunk_track.truth.p.genp + branch_track.truth.p.genp
    else:
        if method == 'range':
            trunk_mom = trunk_track.rangeP.p_muon
            branch_mom = branch_track.rangeP.p_muon
        if method == 'mcs':
            trunk_mom = trunk_track.mcsP.fwdP_muon
            branch_mom = branch_track.mcsP.fwdP_muon
        if method == 'hybrid':
            trunk_mom = trunk_track.mcsP.fwdP_muon
            trunk_mom[TrkInFV(trunk_track.end)] = trunk_track[TrkInFV(trunk_track.end)].rangeP.p_muon
            branch_mom = branch_track.mcsP.fwdP_muon
            branch_mom[TrkInFV(branch_track.end)] = branch_track[TrkInFV(branch_track.end)].rangeP.p_muon
        if method == 'weight_by_len': # Instead of weighting by track momentum which would be most correct, just weight by track length.
            trunk_mom = trunk_track.len
            branch_mom = branch_track.len
        if method == 'dir_only':
            trunk_mom = 1.
            branch_mom = 1.
        #if method == 'track_truth':
        #    trunk_mom = np.sqrt(trunk_track.truth.p.genp.x*trunk_track.truth.p.genp.x +
        #                        trunk_track.truth.p.genp.y*trunk_track.truth.p.genp.y +
        #                        trunk_track.truth.p.genp.z*trunk_track.truth.p.genp.z
        #                       )
        #    branch_mom = np.sqrt(branch_track.truth.p.genp.x*branch_track.truth.p.genp.x +
        #                         branch_track.truth.p.genp.y*branch_track.truth.p.genp.y +
        #                         branch_track.truth.p.genp.z*branch_track.truth.p.genp.z
        #                        )
        scalar_mom = trunk_track.dir.multiply(trunk_mom, axis=0) + \
               branch_track.dir.multiply(branch_mom, axis=0)  

    scalar_dir = scalar_mom.divide(magdf(scalar_mom), axis=0)
    ## 11/27/23: error w/ above line from magdf, but IDK why this worked before. 
    Sbeamangle = np.arccos(scalar_dir.x*beamdir[0] + scalar_dir.y*beamdir[1] + scalar_dir.z*beamdir[2])
    return Sbeamangle

def phi_NuMI(trunk_track, branch_track, beamdir, method): # azimuthal angle around axis from target to ICARUS

    if method == 'p0_p1_truth': #input would be eg evtdf.slc.truth.p0
        scalar_mom = trunk_track.genp + branch_track.genp
    elif method == 'track_truth':
        scalar_mom = trunk_track.truth.p.genp + branch_track.truth.p.genp
    else:
        if method == 'range':
            trunk_mom = trunk_track.rangeP.p_muon
            branch_mom = branch_track.rangeP.p_muon
        if method == 'mcs':
            trunk_mom = trunk_track.mcsP.fwdP_muon
            branch_mom = branch_track.mcsP.fwdP_muon
        if method == 'weight_by_len': # Instead of weighting by track momentum which would be most correct, just weight by track length.
            trunk_mom = trunk_track.len
            branch_mom = branch_track.len
        if method == 'dir_only':
            trunk_mom = 1.
            branch_mom = 1.
        #if method == 'track_truth':
        #    trunk_mom = np.sqrt(trunk_track.truth.p.genp.x*trunk_track.truth.p.genp.x +
        #                        trunk_track.truth.p.genp.y*trunk_track.truth.p.genp.y +
        #                        trunk_track.truth.p.genp.z*trunk_track.truth.p.genp.z
        #                       )
        #    branch_mom = np.sqrt(branch_track.truth.p.genp.x*branch_track.truth.p.genp.x +
        #                         branch_track.truth.p.genp.y*branch_track.truth.p.genp.y +
        #                         branch_track.truth.p.genp.z*branch_track.truth.p.genp.z
        #                        )
        scalar_mom = trunk_track.dir.multiply(trunk_mom, axis=0) + \
               branch_track.dir.multiply(branch_mom, axis=0)    
    
    # Note: I don't know why you get weird distributions for the truth variables - is it real, or is something going wrong with adding the two particles' momenta component-wise?
    
    ## STEP 1: project p into plane perpendicular to BEAMDIR. Call this projection vector pp. (added df since dataframe)
    ppdf = scalar_mom - scalar_mom.dot(BEAMDIR)[:,np.newaxis]*np.repeat([BEAMDIR], repeats=scalar_mom.shape[0], axis=0)
    #if method == 'mcs': print("ppdf.loc[[0]]: ", ppdf.loc[[0]])

    ## STEP 2: Find angle of pp wrt plane made by BEAMDIR and NUMIDIR.
    # Find the normal vector of that plane. Call it N.
    N = np.cross(BEAMDIR, NUMIDIR) # Actually, I will use this as the vector I define phi_NuMI wrt. Then stuff originating in the decay pipe will have phi_NuMI=90deg.

    #phi_NuMI = np.arccos(ppdf.dot(N).divide(np.linalg.norm(N)*magdf(ppdf), axis=0))#*180./math.pi
    
    # use arctan2 method instead:
    # documentation: https://numpy.org/doc/stable/reference/generated/numpy.arctan2.html
    # resulting range is (-180,180)
    x2 = ppdf.dot(N/(np.linalg.norm(N))) # x projection, second arg in arctan2
    x1_proj_axis = np.cross(BEAMDIR,N)
    x1 = ppdf.dot(x1_proj_axis/np.linalg.norm(x1_proj_axis)) # y projection, first arg in arctan2
    phi_NuMI = np.arctan2(x1,x2)
    
    return ppdf, phi_NuMI

def add_calculated_evtdf_cols(df, newcol_name=None, newcol_val=None):
    # Stuff in this function does not have to cross-talk with other dataframes.
    evtdf = df.copy()
    
    # NOTE: Here, "trueTrk" refers to true momentum of trunk/branch TRACK (not neccessarily the muon's)
#       "trueParticle" uses the tracks' truth-matched particle truth information.
    evtdf["Snumi_angle_mcs"] = Sbeamangle(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'mcs')
    evtdf["Snumi_angle_rangeBased"] = Sbeamangle(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'range')
    evtdf["Snumi_angle_hybrid_rangeMCS"] = Sbeamangle(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'hybrid')
    evtdf["Snumi_angle_trkDirOnly"] = Sbeamangle(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'dir_only')
    evtdf["Snumi_angle_planar"] = Sbeamangle(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'planar')
    evtdf["Snumi_angle_wgtByLen"] = Sbeamangle(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'weight_by_len')
    evtdf["Snumi_angle_trueTrk"] = Sbeamangle(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'track_truth')
    evtdf["Snumi_angle_trueParticle"] = Sbeamangle(evtdf.slc.truth.p0, evtdf.slc.truth.p1, BEAMDIR, 'p0_p1_truth')   
    ## Note: NaN for Cosmics b/c p0 and p1 don't mean anything.
    #print(min(evtdf["Snumi_angle_planar"]*180/math.pi), max(evtdf["Snumi_angle_planar"]*180/math.pi), sep=', ')
    evtdf["phi_NuMI_mcs"] = phi_NuMI(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'mcs')[1]
    evtdf["phi_NuMI_rangeBased"] = phi_NuMI(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'range')[1]
    evtdf["phi_NuMI_trueTrk"] = phi_NuMI(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'track_truth')[1]
    evtdf["phi_NuMI_trueParticle"] = phi_NuMI(evtdf.slc.truth.p0, evtdf.slc.truth.p1, BEAMDIR, 'p0_p1_truth')[1]
    
    # add phi_NuMI to evtdf:
    evtdf["phi_NuMI_mcs"] = phi_NuMI(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'mcs')[1]
    evtdf["phi_NuMI_rangeBased"] = phi_NuMI(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'range')[1]
    evtdf["phi_NuMI_trueTrk"] = phi_NuMI(evtdf.trunk.trk, evtdf.branch.trk, BEAMDIR, 'track_truth')[1]
    evtdf["phi_NuMI_trueParticle"] = phi_NuMI(evtdf.slc.truth.p0, evtdf.slc.truth.p1, BEAMDIR, 'p0_p1_truth')[1]
    
    # add opening angle:
    evtdf["opening_angle"] = np.arccos(dotdf(evtdf.trunk.trk.dir, evtdf.branch.trk.dir))
    
    # Add true opening angle to evtdf:
    # WARNING: Take with a grain of salt, as p0,p1 studies in unc_MCstudies nb indicate that truth may be botched for the ALP samples.
    # calculated using the tracks' truth-matched particles' truth information (not the gen-level BSM/parent particle's momentum)
    evtdf["true_opening_angle"] = np.arccos(
        (evtdf.slc.truth.p0.genp.x*evtdf.slc.truth.p1.genp.x +
         evtdf.slc.truth.p0.genp.y*evtdf.slc.truth.p1.genp.y +
         evtdf.slc.truth.p0.genp.z*evtdf.slc.truth.p1.genp.z
        )/(
            np.sqrt(evtdf.slc.truth.p0.genp.x**2 + evtdf.slc.truth.p0.genp.y**2 + evtdf.slc.truth.p0.genp.z**2) *
            np.sqrt(evtdf.slc.truth.p1.genp.x**2 + evtdf.slc.truth.p1.genp.y**2 + evtdf.slc.truth.p1.genp.z**2)
        )
    )
    ## Note: NaN for Cosmics b/c p0 and p1 don't mean anything.
    
    # Add longer/shorter track length to evtdf:
    longer_track_length = []
    shorter_track_length = []
    trunk_len = np.array(evtdf.trunk.trk.len)
    branch_len = np.array(evtdf.branch.trk.len)
    short_is_trunk = 0
    short_is_branch = 0
    for i in range(len(trunk_len)):
        longer_track_length.append(max([ trunk_len[i], branch_len[i] ]))
        shorter_track_length.append(min([ trunk_len[i], branch_len[i] ]))
        if trunk_len[i] < branch_len[i]:
            short_is_trunk = short_is_trunk + 1
        else:
            short_is_branch = short_is_branch + 1
    evtdf["longer_track_length"] = longer_track_length
    evtdf["shorter_track_length"] = shorter_track_length

    # Add columns to dataframe as needed:
    if newcol_name is not None:
        for c, col_name in enumerate(newcol_name):
            evtdf[col_name] = newcol_val[c]
        print("After adding other columns: ", evtdf.shape)
    
    return evtdf

def getp(track, method): # returns a momentum mabnitude.
    if method == 'range': # input would be eg evtdf.trunk.trk
        mom = track.rangeP.p_muon
    if method == 'mcs':
        mom = track.mcsP.fwdP_muon
    if method == 'track_truth':
        mom = np.sqrt(track.truth.p.genp.x*track.truth.p.genp.x +
                      track.truth.p.genp.y*track.truth.p.genp.y +
                      track.truth.p.genp.z*track.truth.p.genp.z
                     )
    if method == 'p0_p1_truth': #input would be eg evtdf.slc.truth.p0
        mom = np.sqrt(track.genp.x*track.genp.x +
                      track.genp.y*track.genp.y +
                      track.genp.z*track.genp.z
                     )
    return mom

# PHYSICS (Constants and branching ratios) - Use this stuff to reweight HPS from kaons to ALPs from kaons. 

# CONSTANTS
# convention: masses in GeV

M_PI = np.pi
elec_mass = 0.0005109989461
muon_mass = 0.1056583745
tau_mass  = 1.77686
piplus_mass = 0.13957039
pizero_mass = 0.1349768
kplus_mass = 0.493677
klong_mass = 0.497611
tquark_mass = 172.76
eta_mass = 0.547862
rho_mass = 0.77526
etap_mass = 0.95778
Wmass = 80.377

fine_structure_constant = 7.2973525693e-3
Gfermi = 1.166379e-5
higgs_vev = 1. / sqrt(sqrt(2)*Gfermi)
sin2thetaW = 0.2312
gL = -0.5 + sin2thetaW
gR = sin2thetaW

fpion = 0.1302 # 130 MeV convention

# compute the eta decay constants
# From https://link.springer.com/article/10.1140/epjc/s10052-021-08861-y
f0 = 0.148
f8 = 0.165
th0 = -6.9*M_PI/180 # rad
th8 = -21.2*M_PI/180. # rad
feta = cos(th8)*f8/sqrt(3) + sin(th0)*f0/sqrt(6) #eq 83
fetap = sin(th8)*f8/sqrt(3) - cos(th0)*f0/sqrt(6) #eq 83

frho = 0.171# GeV^2
grho = 1 - 2*sin2thetaW # // https://link.springer.com/article/10.1140/epjc/s10052-021-08861-y

# unit conversion
hbar = 6.582119569e-16 # GeV*ns or eV*s https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
c_cm_per_ns = 29.9792458 # cm / ns https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf

# kaon lifetimes
kplus_lifetime = 1.238e1 # ns https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
klong_lifetime = 5.116e1 # ns https://pdg.lbl.gov/2020/listings/rpp2020-list-K-zero-L.pdf (FIT)
kshort_lifetime = 8.954e-2 # ns https://pdg.lbl.gov/2014/tables/rpp2014-tab-mesons-strange.pdf

# other lifetimes
tau_lifetime = 290.3e-6 # ns https://pdg.lbl.gov/2019/tables/rpp2019-sum-leptons.pdf

# and widths
rho_width = 0.1478 # GeV https://pdg.lbl.gov/2019/listings/rpp2019-list-rho-770.pdf

# Kaon decay branching ratios
kaonp_mup_numu = 0.6356 # From PDG: https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
kaonp_ep_nue = 1.582e-5 # From PDG: https://pdg.lbl.gov/2020/listings/rpp2020-list-K-plus-minus.pdf
kshort_2pi = 0.692 # https://pdg.lbl.gov/2014/tables/rpp2014-tab-mesons-strange.pdf

# CKM matrix
abs_Vud_squared = 0.97370 * 0.97370 # https://pdg.lbl.gov/2020/reviews/rpp2020-rev-ckm-matrix.pdf (12.7)

abs_VtsVtd_squared = 1.19777e-7
rel_VtsVtd_squared = 1.02136e-7

def child_momentum(M, m1, m2) :
    return sqrt(M * M * M * M
        -2 * M * M * m1 * m1
        -2 * M * M * m2 * m2
           + m1 * m1 * m1 * m1
           + m2 * m2 * m2 * m2
        -2 * m1 * m1 * m2 * m2) / ( 2 * M )

higgs_momentum = child_momentum

# SCALAR WIDTHS

# Kpm -> S
def HPS_KaonPlusBranchingRatio(higs_mass, mixing):
    # kplus width
    #
    # matrix eement for kplus
    M_KP = (1. / 2.) * ( 3. / (16. * M_PI * M_PI * higgs_vev * higgs_vev * higgs_vev)) * (kplus_mass * kplus_mass) * (tquark_mass * tquark_mass) * mixing
    M_KP2 = M_KP * M_KP * abs_VtsVtd_squared

    kplus_width = (2 * higgs_momentum(kplus_mass, piplus_mass, higs_mass)/kplus_mass) * M_KP2 / (16 * M_PI * kplus_mass)

    # convert to partia ifetime
    kplus_2higgs_lifetime = hbar / kplus_width 

    # and branching ratio
    #
    # (this higgs decay woud make a negigibe contribution to the overa ifetime)
    return kplus_lifetime / kplus_2higgs_lifetime
        

# KL -> S
def HPS_KaonLongBranchingRatio(higs_mass, mixing):
    M_KL = (1. / 2.) * (3. / (16. * M_PI * M_PI * higgs_vev * higgs_vev * higgs_vev)) * (klong_mass * klong_mass) * (tquark_mass * tquark_mass) * mixing
    M_KL2 = M_KL * M_KL * rel_VtsVtd_squared

    klong_width = (2 * higgs_momentum(klong_mass, pizero_mass, higs_mass) / klong_mass) * M_KL2 / (16 * M_PI * klong_mass)

    klong_2higgs_lifetime = hbar / klong_width

    return klong_lifetime / klong_2higgs_lifetime

# Get the partial width for lepton decays
def LeptonPartialWidth(lep_mass, higs_mass, mixing):
    if (lep_mass * 2 >= higs_mass):
        width = 0.
    else:
        width = (mixing * mixing * lep_mass * lep_mass * higs_mass / (8 * M_PI * higgs_vev * higgs_vev)) * pow(1 - 4 * lep_mass * lep_mass / (higs_mass * higs_mass), 3. / 2.)
    #width[(lep_mass * 2 >= higs_mass)] = 0
    return width

# S -> ee
def HPS_ElectronPartialWidth(higs_mass, mixing):
    return LeptonPartialWidth(elec_mass, higs_mass, mixing)

# S -> mumu
def HPS_MuonPartialWidth(higs_mass, mixing):
    return LeptonPartialWidth(muon_mass, higs_mass, mixing)

# S -> pi pi
def HPS_PionPartialWidth(pion_mass, higs_mass, mixing):
    form_factor = (2. / 9.) * higs_mass * higs_mass + (11. / 9.) * pion_mass * pion_mass
    if (pion_mass * 2 >= higs_mass):
        width = 0.
    else:
        width = (mixing * mixing * 3 * form_factor * form_factor / (32 * M_PI * higgs_vev * higgs_vev * higs_mass)) * pow(1- 4. * pion_mass * pion_mass / (higs_mass * higs_mass), 1. / 2.)
    #width[(pion_mass * 2 >= higs_mass)] = 0
    return width


def HPS_PiPlusPartialWidth(higs_mass, mixing):
    return HPS_PionPartialWidth(piplus_mass, higs_mass, mixing)


def HPS_PiZeroPartialWidth(higs_mass, mixing):
    return HPS_PionPartialWidth(pizero_mass, higs_mass, mixing)

def HPS_total_decay_Width(higs_mass, mixing):
    return(
        HPS_ElectronPartialWidth(higs_mass, mixing) +
        HPS_MuonPartialWidth(higs_mass, mixing) +
        HPS_PiPlusPartialWidth(higs_mass, mixing) + 
        HPS_PiZeroPartialWidth(higs_mass, mixing)
    )

# AXION WIDTHS

# Axion strong production
def ALP_KPlusStrongBranchingRatio(mass, fa, c3):
    ks_partial_lifetime = kshort_lifetime / kshort_2pi
    mass_factor = ((kplus_mass**2 - mass**2)/(4*kplus_mass**2 - 3*mass**2 - piplus_mass**2))**2
    momentum_factor = child_momentum(kplus_mass, mass, piplus_mass) / sqrt(kplus_mass**2 - 4*piplus_mass**2)
    
    lifetime = ks_partial_lifetime / ((4*fpion**2*c3**2/fa**2)*mass_factor*momentum_factor)
    
    if mass > (kplus_mass - piplus_mass):
        BR = 0.
    else:
        BR = kplus_lifetime/lifetime
    #BR[mass > kplus_mass - piplus_mass] = 0
    return BR

def ALP_KLongStrongBranchingRatio(mass, fa, c3): # same as ALP_KPlusStrongBranchingRatio (see Slack conversation with Josh 2/20/24)
    ks_partial_lifetime = kshort_lifetime / kshort_2pi
    mass_factor = ((kplus_mass**2 - mass**2)/(4*kplus_mass**2 - 3*mass**2 - piplus_mass**2))**2
    momentum_factor = child_momentum(kplus_mass, mass, piplus_mass) / sqrt(kplus_mass**2 - 4*piplus_mass**2)
    
    lifetime = ks_partial_lifetime / ((4*fpion**2*c3**2/fa**2)*mass_factor*momentum_factor)
    
    if (mass > kplus_mass - piplus_mass):
        BR = 0.
    else:
        BR = klong_lifetime/lifetime
    #BR[mass > kplus_mass - piplus_mass] = 0
    return BR
    
    
def ALP_KPlusWeakBranchingRatio(mass, fa, c2):
    def f(x):
        return x*(1+x*(np.log(x) - 1)) / (1-x)**2
    
    # weak coupling constant
    aW = fine_structure_constant / sin2thetaW
    
    ckm_factor = abs_VtsVtd_squared*f(tquark_mass**2/Wmass**2)**2
    
    gasd2 = ckm_factor*(3*sqrt(2)*Gfermi*Wmass**2*c2*aW/(32*M_PI**3*fa))**2
    
    print(aW, gasd2)
    
    width = ((kplus_mass**2)/(32*M_PI))*(1-piplus_mass**2/kplus_mass**2)**2*gasd2*child_momentum(kplus_mass, mass, piplus_mass)
    
    lifetime = hbar / width
    
    return kplus_lifetime / lifetime

def ALP_KLongWeakBranchingRatio(mass, fa, c2):
    def f(x):
        return x*(1+x*(np.log(x) - 1)) / (1-x)**2
    
    # weak coupling constant
    aW = fine_structure_constant / sin2thetaW
    
    ckm_factor = (abs_VtsVtd_squared - rel_VtsVtd_squared)*f(tquark_mass**2/Wmass**2)**2
    
    gasd2 = ckm_factor*(3*sqrt(2)*Gfermi*Wmass**2*c2*aW/(32*M_PI**3*fa))**2
        
    width = ((klong_mass**2)/(32*M_PI))*(1-pizero_mass**2/klong_mass**2)**2*gasd2*child_momentum(klong_mass, mass, pizero_mass)
    
    lifetime = hbar / width
    
    return klong_lifetime / lifetime

def ALP_uu_width(mass, fa, cAu):
    return float( (cAu**2*mass*muon_mass**2 / (8*np.pi*fa**2))*sqrt(1-4*muon_mass**2/mass**2) )

def ALP_gg_width(mass, fa, c1, c2, c3, cAu):
    def f(x):
        xl1 = np.pi/2 - (1j/2)*np.log(1+np.sqrt(1-x)/(1-np.sqrt(1-x)))
        if x < 1:
            return xl1
        else:
            xg1 = np.arcsin(1/np.sqrt(x))
            xret = xg1.copy()
            return xret
        #xret[x < 1] = xl1[x < 1] # TODO: I (Jamie) commented this out. Do I need it?
        
    def B(x):
        return 1-x*f(x)**2
    
    
    cg = B(4*muon_mass**2/mass**2)*cAu + (5/3)*c1 + c2 + c3*(-1.93 + (1/3)*mass**2/(mass**2 - pizero_mass**2) +\
                     (8/9)*(mass**2 - 4*pizero_mass**2/9)/(mass**2 - eta_mass**2) +\
                     (7/9)*(mass**2-16*pizero_mass**2/9)/(mass**2 - etap_mass**2))
    cg2 = np.abs(cg)**2
    
    return fine_structure_constant**2*cg2*mass**3/(256*M_PI**3*fa**2)

def ALP_total_decay_Width(mass, fa, cAu, c1, c2, c3): # JD: Does not include hadronic decays. Decay to pions matters above about 400MeV mass for the ALPs. TODO (if I ever use this): include anything in the alp_make_decay tool in MeVPrtl
    if mass < 3*pizero_mass:
        return( ALP_uu_width(mass, fa, cAu) + ALP_gg_width(mass, fa, c1, c2, c3, cAu) )
    else:
        print("You can't use this function to get total width! You need to also consider hadronic decays.")
        return -1
    
    
