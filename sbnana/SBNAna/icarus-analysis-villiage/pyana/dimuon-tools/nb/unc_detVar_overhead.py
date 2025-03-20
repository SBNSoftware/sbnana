# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~ Stuff for detector variation studies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In unc_MC_overhead, I make one evtdf for all the nominal samples.
# Here, make an evtdf for each detector variation sample.

import pandas as pd
import pickle
import math
from unc_MC_overhead import *

sets_of_detVar_files_by_sample = [
    higgs_220_detVar_files,
    higgs_240_detVar_files,
    higgs_260_detVar_files,
    higgs_280_detVar_files,
    higgs_300_detVar_files,
    higgs_340_detVar_files,
    cohlike_nu_detVar_files
]
labels_for_samples_w_dvs = labels[:6]+[labels[-1]]
detVar_evtdf = [] # one evtdf for each type of detector variation.
for dv, dvl in enumerate(detVar_labels):
    temp_files = [myset[dv] for myset in sets_of_detVar_files_by_sample]
    temp_mchdfs = [pd.read_hdf(f, key="mch") for f in temp_files]
    temp_mcdfs = [pd.read_hdf(f, key="mcnu") for f in temp_files]
    temp_masses = [int(round(df.iloc[(0)].M*1000.)) for df in temp_mchdfs[:-1]]+[-1]
    th_or_fa_or_null = [float(df.iloc[(0)].C1) for df in temp_mchdfs[:-1]]+[-1]
    cAl_or_null = [float(df.iloc[(0)].C1) for df in temp_mchdfs[:-1]]+[-1]

    temp_evtdfs = [pd.read_hdf(f, key="evt") for f in temp_files]
    temp_hdrs = [pd.read_hdf(f, key="hdr") for f in temp_files]
    temp_pots_notPreferred = [np.sum(hdr.pot * hdr.first_in_subrun) for hdr in temp_hdrs]
    temp_pots = temp_pots_notPreferred

    for i, df in enumerate(temp_evtdfs):
        df["scale"] = GOAL_POT / temp_pots[i]
        df["bsm_mass"] = temp_masses[i]
        df["sample"] = labels_for_samples_w_dvs[i]
        df["higgs"] = False
        df["alp"] = False
        df["nu"] = False
        if '_S' in labels_for_samples_w_dvs[i]: df["higgs"] = True
        if 'alp' in labels_for_samples_w_dvs[i]: df["alp"] = True
        if 'nu' in labels_for_samples_w_dvs[i]: df["nu"] = True
    temp_evtdf = jamie_sample_concat(temp_evtdfs)
    print('evtdf.shape before enforcing correct FV w/ uncontained track lenght:', temp_evtdf.shape)
    temp_evtdf = temp_evtdf[satisfies_new_FV(temp_evtdf)]
    temp_evtdf = add_calculated_evtdf_cols(temp_evtdf)
    detVar_evtdf.append(temp_evtdf)



