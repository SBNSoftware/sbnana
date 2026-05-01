# Useful for:
#    unc_MC Event Selection.ipynb
#    unc_MCstudies.ipynb

# Here, I use the pandas function ' pd.read_hdf(filename, key='evt')' to get the evt dataframes from the files listed in unc_samples.py. This keeps the dataframes fairly light (i.e. they do not have all the weights needed to run systematics on them), so I can handle the full dataset. 

import pandas as pd
import pickle
import math
import var
import re
from data import *
from unc_samples import *

# ~~~~~~~~~~ Stuff from unc_MC Event Selection NB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HPS

higgs_mcdfs = [pd.read_hdf(f, key="mch") for f in higgs_files]
higgs_mcnudfs = [pd.read_hdf(f, key="mcnu") for f in higgs_files]
for i, mch in enumerate(higgs_mcdfs):
    Sdir = mch.exit-mch.enter
    Sdir_mag = np.linalg.norm(Sdir, axis=1)
    Sdir = Sdir/Sdir_mag[:, None]
    mch["POIdir","x"] = Sdir.x # "POI" is for "particle of interest" (scalar, alp, or nu)
    mch["POIdir","y"] = Sdir.y
    mch["POIdir","z"] = Sdir.z
    higgs_mcnudfs[i]["POIdir_x"] = Sdir.x
    higgs_mcnudfs[i]["POIdir_y"] = Sdir.y
    higgs_mcnudfs[i]["POIdir_z"] = Sdir.z
#for c in higgs_mcdfs[0].columns: print(c)
higgs_masses = [int(round(df.iloc[(0)].M*1000.)) for df in higgs_mcdfs]
higgs_thetas = [float(df.iloc[(0)].C1) for df in higgs_mcdfs]
higgs_mass_labels = []
higgs_th_labels = []
higgs_labels = []
for i in range(len(higgs_files)):
    higgs_mass_labels.append("$M_S$ = %a" % higgs_masses[i])
    higgs_th_labels.append("$\\theta_S$ = %a" % higgs_thetas[i])
    higgs_labels.append(higgs_mass_labels[i] + ", " + higgs_th_labels[i])

# ALP Files (ALPs WITHOUT MASS SUPPRESSION)

alp_nosup_mcdfs = [pd.read_hdf(f, key="mch") for f in alp_nosup_files]
alp_nosup_mcnudfs = [pd.read_hdf(f, key="mcnu") for f in alp_nosup_files]
for i, mch in enumerate(alp_nosup_mcdfs):
    ALPdir = mch.exit-mch.enter
    ALPdir_mag = np.linalg.norm(ALPdir, axis=1)
    ALPdir = ALPdir/ALPdir_mag[:, None]
    mch["POIdir","x"] = ALPdir.x
    mch["POIdir","y"] = ALPdir.y
    mch["POIdir","z"] = ALPdir.z
    alp_nosup_mcnudfs[i]["POIdir_x"] = ALPdir.x
    alp_nosup_mcnudfs[i]["POIdir_y"] = ALPdir.y
    alp_nosup_mcnudfs[i]["POIdir_z"] = ALPdir.z

alp_nosup_masses = [int(round(df.iloc[(0)].M*1000.)) for df in alp_nosup_mcdfs]
alp_nosup_fa = [float(df.iloc[(0)].C1) for df in alp_nosup_mcdfs] #play with this one!
alp_nosup_cAl = [float(df.iloc[(0)].C2) for df in alp_nosup_mcdfs]
#print(alp_nosup_masses, alp_nosup_fa, alp_nosup_cAl, sep='\n')
alp_nosup_mass_labels = []
alp_nosup_fa_labels = []
alp_nosup_cAl_labels = []
alp_nosup_labels = []
for i in range(len(alp_nosup_files)):
    alp_nosup_mass_labels.append("$M_{ALP}$ = %a" % alp_nosup_masses[i])
    alp_nosup_fa_labels.append("$fa$ = %a" % "{:.1e}".format(alp_nosup_fa[i]))
    if len(str(alp_nosup_fa[i]))-str(alp_nosup_fa[i]).count('0') > 2:
        print('WARNING: Your fa labels might be misleading. You should double check them, and fix if needed.')
    alp_nosup_cAl_labels.append("$c$ = %a" % alp_nosup_cAl[i])
    alp_nosup_labels.append(alp_nosup_mass_labels[i] + ", " + alp_nosup_cAl_labels[i] + ", " + alp_nosup_fa_labels[i])

# NEUTRINOS

#mcfile = "/icarus/data/users/gputnam/DMCP2023G/mc/dfs-wstubs/mcnuphase2_evt_reprodD.df"

#nu_files = [nu_file, cohlike_nu_file, cohlike_nu_file_extra_2501, 
#            cohlike_nu_file_extra_2502_12, cohlike_nu_file_extra_2502_34,
#            cohlike_nu_file_extra_2502_5, 
#            cohlike_nu_file_extra_2502_6789]
#nu_labels = ["$\\nu$", "Coh-like $\\nu$", "Coh-like $\\nu$ (ad'l 2501)", 
#             "Coh-like $\\nu$ (ad'l 2502_12)", "Coh-like $\\nu$ (ad'l 2502_34)", #"Coh-like $\\nu$ (ad'l 2502_5)",
#             "Coh-like $\\nu$ (ad'l 2502_6789)"]#*len(nu_files)

#nu_files = [nu_file,
#            cohlike_nu_file_extra_2501, 
#            cohlike_nu_file_extra_2502_12, 
#            cohlike_nu_file_extra_2502_5, 
#            cohlike_nu_file_extra_2502_6789]
#nu_labels = ["$\\nu$",
#            "Coh-like $\\nu$ (ad'l 2501)", 
#             "Coh-like $\\nu$ (ad'l 2502_12)",
#             "Coh-like $\\nu$ (ad'l 2502_5)",
#             "Coh-like $\\nu$ (ad'l 2502_6789)"]#*len(nu_files)

nu_files = [nu_file,
            cohlike_1,
            cohlike_2
            ]
nu_labels = ["$\\nu$",
             "Coh-like $\\nu$ (F1)", 
             "Coh-like $\\nu$ (F2)"
             ]

nu_masses = [-1]*len(nu_files)
nu_mcnudfs = [pd.read_hdf(f, key="mcnu") for f in nu_files]
for df in nu_mcnudfs:
    nudir = df.momentum
    nudir_mag = np.linalg.norm(nudir, axis=1)
    nudir = nudir/nudir_mag[:, None]
    df["POIdir_x"] = nudir.x
    df["POIdir_y"] = nudir.y
    df["POIdir_z"] = nudir.z

# PUT IT ALL TOGETHER:

labels = higgs_labels + alp_nosup_labels + nu_labels
masses = higgs_masses + alp_nosup_masses + nu_masses # MeV
df_files = higgs_files + alp_nosup_files + nu_files
mcnudfs = higgs_mcnudfs + alp_nosup_mcnudfs + nu_mcnudfs

# REFORMAT THE LABELS:
test_labels = []
for i, l in enumerate(labels):
    print(l)
    new_label = l
    new_label = new_label.replace("$c$","$c_{\mu}$")
    new_label = new_label.replace("$fa$","$f_a$")
    new_label = new_label.replace("'","")
    new_label = new_label.replace("1.0","1")
    new_label = new_label.replace("'","")
    
    new_label_split = new_label.split(',')
    print(new_label_split)
    for s, sp in enumerate(new_label_split):
        if 'M_' in sp:
            new_label_split[s] = sp+' MeV'
        if 'f_a' in sp:
            new_label_split[s] = sp+' GeV'
    new_label=", ".join(new_label_split)
    
    match = re.search("e[\+\-]\d+", new_label) 
    if match:
        etothe = match.group()
        exponent = re.search("0\d+", etothe)
        if exponent:
            exponent = exponent.group()[1:]
        new_scinot = re.sub("0\d+", exponent, etothe, count=0, flags=0)
        new_scinot = 'e{'+new_scinot[1:]+'}'
        new_scinot = new_scinot.replace('e','*10^')
        new_scinot = new_scinot.replace('+','')
        new_scinot = "$"+new_scinot+"$"
        #print(new_scinot)
        #print(new_label, ', ', etothe, ', ', new_scinot)
        new_label = new_label.replace(etothe,new_scinot)
        new_label = new_label.replace('1$*10','$10')
        new_label = new_label.replace('*','\\times')
    test_labels.append(new_label)
labels = test_labels

## IF YOU DO NOT CARE ABOUT CV WEIGHTS:

#evtdfs = [pd.read_hdf(f, key="evt") for f in df_files]

## IF YOU DO CARE ABOUT CV WEIGHTS (***): #takes 30 minutes.

def add_cv_weights(file, mccut=None,
                   mc=True, mc_incoh=False, mccoh=False,
                   bsm=False, hps=False, alp=False,
                   onbeam=False, offbeam=False, Run1=False, Run2=False,
                  ):
    evt = mc_dataset(file, "evt", mcnukey="mcnuwgt", syst_weights=True, isscalar=bsm, alp=alp)
    #evt = mc_dataset(file, "evt", mcnukey="mcnuwgt", syst_weights=False)#, mccut=mccut)
    evt_df = evt.df
    if mccut is not None:
        print('recognized mccut')
        print('shape of evt_df before applying mccut: ', evt_df.shape)
        evt_df = evt_df[mccut(evt_df)]
        print('shape of evt_df after applying mccut: ', evt_df.shape)
    evt_df["mc"] = mc
    evt_df["mc_incoh"] = mc_incoh
    evt_df["mccoh"] = mccoh
    evt_df["bsm"] = bsm
    evt_df["hps"] = hps
    evt_df["alp"] = alp
    evt_df["onbeam"] = evt_df["offbeam"] = evt_df["Run1"] = evt_df["Run2"] = False
    evt_df['livetime'] = evt.livetime
    evt_df['POT'] = evt.POT
    #evt_df['scale'] = GOAL_POT / evt.POT
    #print('livetime: ', evt.livetime)
    #print('POT: ', evt.POT)
    
    evt_df = evt_df[satisfies_new_FV(evt_df)] # Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events). This requires that any uncontained tracks are at least 1m long.
    when_uncontained = ~TrkInFV(evt_df.trunk.trk.end) | ~TrkInFV(evt_df.branch.trk.end) # Require exiting. (Min track length of 1m for exiting track already accounted for in satisfies_new_FV in line above.)
    evt_df = evt_df[when_uncontained]
    evt_df = add_calculated_evtdf_cols(evt_df)
    
    return evt_df

higgs_evtdfs = [add_cv_weights(file, bsm=True, hps=True) for file in higgs_files]
alp_evtdfs = [add_cv_weights(file, bsm=True, alp=True) for file in alp_nosup_files]
nu_evtdfs = [add_cv_weights(nu_files[0], mccut=is_NOT_coh_like_JD, mc_incoh=True),
             add_cv_weights(nu_files[1], mccut=is_coh_like_JD, mccoh=True),
             add_cv_weights(nu_files[2], mccut=is_coh_like_JD, mccoh=True)
            ]
evtdfs = higgs_evtdfs + alp_evtdfs + nu_evtdfs
evtdf = jamie_sample_concat(evtdfs)

## list of files that have cv weights in them (but not all the other heavy stuff for systematics)
#for f in df_files:
#    print(f)
#    print(f.split('/'))
#w_cv_weights_dir = '/exp/icarus/data/users/jdyer/muon_selection/bsm_sample_dfs_w_cv_weights/'
#bsm_cv_df_files = [w_cv_weights_dir+f.split('/')[-1] for f in df_files]
#bsm_cv_df_files = bsm_cv_df_files[:-2]
#print('')
#for f in bsm_cv_df_files: print(f)

## DO THIS ONCE AND SAVE: needed to include cv weight for concrete and Minerva tune
##evtdfs = []
#pots = []
#for f in bsm_cv_df_files:
#    print('\n \n \n \n ', f, '\n')
#    dataset = data.mc_dataset(f, "evt", mcnukey='mcnuwgt')
#    pots.append(dataset.POT)
#    temp_df = dataset.df
#    cols = list(temp_df.columns)
#    cv_ind = cols.index(('wgt', 'cv', '', '', '', ''))
#    assert cols[cv_ind+1] == ('wgt', 'concrete', 'cv', '', '', '')
#    assert cols[cv_ind+2] == ('wgt', 'coh', 'cv', '', '', '')
#    temp_df = temp_df.drop(columns=cols[(cv_ind+3):]) # delete the unneeded columns to lighten the load
#    #evtdfs.append(temp_df)
#    with pd.HDFStore(f) as hdf:
#        hdf.put(key='evt_w_cv', value=temp_df)
#    del temp_df
##    # 24/05/16 NOTE: My alp events are not working with data.mc_dataset function, which I need to get cv weights. So for now just (misleadingly) copy a Higgs dataframe there.
#for (i,f) in enumerate(bsm_cv_df_files):
#    print(f, ': ', pots[i], ' POT')

## OTHERWISE (IF YOU'VE ALREADY DONE THE ABOVE THING ONCE):
# Note: this is a bit Frankensteined (and unfortunately takes a while to run the
#       data.mc_dataset function on the neutrinos) b/c not enough space in my data area to copy/resave the nu dataframes
#       with the cv weights added in (and don't have write access to the ones in Gray's area.)
#evtdfs = [pd.read_hdf(f, key="evt_w_cv") for f in bsm_cv_df_files]

#for nu_file in nu_files:
#    dataset = data.mc_dataset(nu_file, "evt", mcnukey='mcnuwgt', syst_weights=False)
#    temp_df = dataset.df
#    cols = list(temp_df.columns)
#    cv_ind = cols.index(('wgt', 'cv', '', '', '', ''))
#    assert cols[cv_ind+1] == ('wgt', 'concrete', 'cv', '', '', '')
#    assert cols[cv_ind+2] == ('wgt', 'coh', 'cv', '', '', '')
#    temp_df = temp_df.drop(columns=cols[(cv_ind+3):])
#    print(nu_file, ': ', dataset.POT, ' POT')
#    evtdfs.append(temp_df)
#    del temp_df

# END OF ***

## MASK OUT COHLIKE EVENTS IN THE NOMINAL NU SAMPLE:
## Made this fix 2/8/25, updated 2/19/25
## NOTE 7/29/25: commented out today because I just remade dataframes with CAFs on hand, and don't currently have nominal nu dataframe.
#nomnu = evtdfs[-7]
#new_n_incoh = nomnu[~is_coh_like_JD(nomnu)].shape[0]
#print('Percentage of non-cohlike events in the nominal nu sample: ', (new_n_incoh)/nomnu.shape[0])
#evtdfs[-7] = nomnu[~is_coh_like_JD(nomnu)]


hdrs = [pd.read_hdf(f, key="hdr") for f in df_files]
pots_notPreferred = [np.sum(hdr.pot * hdr.first_in_subrun) for hdr in hdrs]
pots = pots_notPreferred # Note: May want to double check that this is consistent with text block below (but I expect it to always be, for MC).
#mcnudfs = [pd.read_hdf(f, key="mcnu") for f in df_files]

# Make sure your POT-accounting is correct for cohlike events since not all generated at once as a single sample

#total_cohlike_pot = pots[-6] + pots[-5] + pots[-4] + pots[-3] + pots[-2] + pots[-1]
#total_cohlike_pot = pots[-4] + pots[-3] + pots[-2] + pots[-1]
total_cohlike_pot = pots[-2] + pots[-1]
#pots[-6] = total_cohlike_pot
#pots[-5] = total_cohlike_pot
#pots[-4] = total_cohlike_pot
#pots[-3] = total_cohlike_pot
pots[-2] = total_cohlike_pot
pots[-1] = total_cohlike_pot

for i in range(len(df_files)):
    evtdfs[i]["scale"] = GOAL_POT / pots[i]
    #mcdfs[i]["scale"] = GOAL_POT / pots[i]

    evtdfs[i]["bsm_mass"] = masses[i]
    evtdfs[i]["sample"] = labels[i]
#    evtdfs[i]["higgs"] = False
#    evtdfs[i]["alp"] = False
#    evtdfs[i]["nu"] = False
#    evtdfs[i]["mc_incoh"] = False
#    evtdfs[i]["mccoh"] = False
#    x = df_files[i].split('/')
#    #print(x)
#    if ("Higgs" in x[-1]) | ("ms" in x[-1]):
#        evtdfs[i]["higgs"] = True
#    if ("alp" in x[-1]) | ("ma" in x[-1]):
#        evtdfs[i]["alp"] = True
#    if ("Nu" in x[-1]) | ("Coh" in x[-1]):
#        evtdfs[i]["nu"] = True
#        if ("Nu" in x[-1]):
#            evtdfs[i]["mc_incoh"] = True
#        if ("Coh" in x[-1]):
#            evtdfs[i]["mccoh"] = True
            
evtdf = jamie_sample_concat(evtdfs)
print('evtdf.shape before enforcing correct FV w/ uncontained track length:', evtdf.shape)
evtdf = evtdf[satisfies_new_FV(evtdf)]
print('evtdf.shape after enforcing correct FV w/ uncontained track length:', evtdf.shape)

# ADD SOME STUFF TO EVTDF:

# ~~ stuff below replaced with a function September 12, 2024 ~

evtdf = add_calculated_evtdf_cols(evtdf) # THIS MUST BE CALLED AFTER REQUIRING EXITING FOR THE NDAUGHTERS FOR CONTAINED TRACKS VARIABLE TO WORK!!!

# Here, I have made one evtdf for all the nominal samples. 
# I make dataframes for the detector variation samples in unc_detVar_overhead.py.


