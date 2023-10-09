# functions copied from unc_MC Event Selection.ipynb on September 28, 2023

from math import dist
import numpy as np

def decay_weight(in_dist, out_dist, mean_dist):
    return np.exp(-in_dist / mean_dist) - \
          np.exp(-out_dist / mean_dist)
    
def flux_weight(mixing):
    return mixing * mixing

def reweight_mixing(newmixing, start, enter, exit, mean_dist, oldmixing=1e-5):
    dist_in = dist(enter, start)
    dist_out = dist(exit, start)
    
    old_dcy = decay_weight(dist_in, dist_out, mean_dist)
    new_dcy = decay_weight(dist_in, dist_out, mean_dist / (newmixing**2 / oldmixing**2))
    
    old_flux = flux_weight(oldmixing)
    new_flux = flux_weight(newmixing)
    
    return (new_flux / old_flux) * (new_dcy / old_dcy)



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
    
def add_hdr_info(df):
    new_df = df.copy()
    proc = []
    clust = []
    for idx in final_df.index:
        p = hdrs[idx[0]].loc[(idx[1],idx[2])].proc
        c = hdrs[idx[0]].loc[(idx[1],idx[2])].cluster
        proc.append(p)
        clust.append(c)
    new_df['proc'] = proc
    new_df['cluster'] = clust
    return new_df

# COLORS

blues = ["#B0E0E6", "#87CEEB", "#6495ED", "#1E90FF","#4682B4", "#00008B"]
greens = ["#CFFCDA", "#8DF9A7", "#72DF8C", "#48B161", "#2C8942", "#145C25", "#054915"]
purples = ["#D8BFD8", "#9370DB", "#663399"]
oranges = ["#FECC88", "#FCC171", "#ECA646", "#DE9025", "#CD7A08", "#A15F03"]#, "#", "#"]
# Use this to find html hex string color codes: https://htmlcolorcodes.com 

# BROAD CATEGORIES ( BSM model + benchmark | nu | cosmic )
# tmatch = truth match (based on what contributed most E to slice)

def make_categories(df, bsm=False): 
    
    is_higgs = (df.slc.tmatch.idx >= 0) & df.higgs & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0) # Only consider muon channel, exclude any Higgs that decayed to pions.
    #is_higgs.name = "Scalar"
    higgs_benchmarks = []
    higgs_benchmark_names = df[is_higgs]["sample"].unique()
    for l in range(len(higgs_benchmark_names)):
        #print(l)
        bm = is_higgs & (df["sample"] == higgs_benchmark_names[l])
        bm.name = higgs_benchmark_names[l]
        bm.color = blues[l]
        higgs_benchmarks.append( bm )    
    
    is_alp = (df.slc.tmatch.idx >= 0) & df.alp_withsup & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0)
    #is_alp.name = "ALP"
    alp_benchmarks = []
    alp_benchmark_names = df[is_alp]["sample"].unique()
    for l in range(len(alp_benchmark_names)):
        bm = is_alp & (df["sample"] == alp_benchmark_names[l])
        bm.name = alp_benchmark_names[l]
        bm.color = purples[l]
        alp_benchmarks.append( bm )
        
    is_alp_nosup = (df.slc.tmatch.idx >= 0) & df.alp_nosup & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0)
    #is_alp.name = "ALP"
    alp_nosup_benchmarks = []
    alp_nosup_benchmark_names = df[is_alp_nosup]["sample"].unique()
    for l in range(len(alp_nosup_benchmark_names)):
        bm = is_alp_nosup & (df["sample"] == alp_nosup_benchmark_names[l])
        bm.name = alp_nosup_benchmark_names[l]# + " (no sup.)"
        bm.color = greens[l]
        alp_nosup_benchmarks.append( bm )
        
    #is_nu = (df.slc.tmatch.idx >= 0) & ~df.higgs
    is_nu = (df.slc.tmatch.idx >= 0) & df.nu
    is_nu.name = "$\\nu$"
    is_nu.color = "C1"
    
    is_cosmic = df.slc.tmatch.idx < 0 # from any file, so I'm assuming both samples include cosmics.
    is_cosmic.name = "Cosmic"
    is_cosmic.color = "C3"

    if bsm:
        is_bsm = (df.slc.tmatch.idx >= 0) & (df.slc.truth.npi == 0) & (df.slc.truth.npi0 == 0) & ( df.higgs | df.alp_withsup | df.alp_nosup )
        is_bsm.name = "BSM"
        is_bsm.color = "#191970" #, "#00FF7F" #"#000000"
        categories = [is_bsm] + [is_nu] + [is_cosmic]
        return categories
    
    else:
        #categories = [is_higgs, is_axion, is_alp, is_nu, is_cosmic]
        categories = higgs_benchmarks + alp_benchmarks + alp_nosup_benchmarks + [is_nu] + [is_cosmic]
        return categories

#categories = make_categories(evtdf)

