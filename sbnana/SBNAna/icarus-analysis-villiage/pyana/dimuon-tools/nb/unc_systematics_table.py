import numpy as np
import pandas as pd

output_path = '/exp/icarus/data/users/jberger/dimuon_analysis_031026/systematics_output/'
selected_filename = 'selected_events.pkl'


try:
    selected_evtdf = pd.read_pickle(output_path + selected_filename)
except FileNotFoundError:
    from unc_MC_overhead import evtdf,hdrs
    from unc_funcs import add_hdr_info
    from unc_cuts import evtSel_dict,apply_cuts

    cut_list = evtSel_dict['DD'][0]
    thresholds = evtSel_dict['DD'][1]
    cut_results = apply_cuts(evtdf, cut_list, thresholds = thresholds, detailed_nu='none', detailed_bsm=True)
    mask = cut_results[-1]
    selected_evtdf = add_hdr_info(evtdf[mask], hdrs)
    selected_evtdf.to_pickle(output_path + selected_filename)

def print_df(df,name):
    tot_evt = (df.scale *df.wgt.cv.tot).sum()
    df_trkrwt = np.maximum(df[('wgt','cv','cathode_crossing','','','')] * df[('wgt','cv','z_crossing','','','')],0)
    univ_tot = [(df.scale * df.wgt.cv.tot * np.divide(df[("wgt","track_split","univ_" + str(i),"","","")],df_trkrwt)).sum() for i in range(100)]
    tot_err = np.std(univ_tot)
    tot_rel = tot_err / tot_evt
    row = [name, len(df), np.round(tot_evt,3), np.round(tot_err,3), np.round(tot_rel,3)]
    print('{: <65} {: >8} {: >8} {: >8} {: >8}'.format(*row))


is_higgs = ( (selected_evtdf.slc.tmatch.idx >= 0) & selected_evtdf.hps & (selected_evtdf.slc.truth.npi == 0) & (selected_evtdf.slc.truth.npi0 == 0) )# Only consider muon channel, exclude any Higgs that decayed to pions. #250919
higgs_benchmark_names = selected_evtdf[is_higgs]["sample"].unique()
for bench in higgs_benchmark_names:
    bench_df = selected_evtdf[is_higgs][selected_evtdf[is_higgs]["sample"] == bench]
    print_df(bench_df,bench)

is_alp = ( (selected_evtdf.slc.tmatch.idx >= 0) & selected_evtdf.alp & (selected_evtdf.slc.truth.npi == 0) & (selected_evtdf.slc.truth.npi0 == 0) )# Only consider muon channel, exclude any Higgs that decayed to pions. #250919
alp_benchmark_names = selected_evtdf[is_alp]["sample"].unique()
for bench in alp_benchmark_names:
    bench_df = selected_evtdf[is_alp][selected_evtdf[is_alp]["sample"] == bench]
    print_df(bench_df,bench)

is_incoh = (selected_evtdf.slc.tmatch.idx >= 0) & (selected_evtdf.mc_incoh)
bench = "nu except coherent"
bench_df = selected_evtdf[is_incoh]
print_df(bench_df,bench)

is_coh = (selected_evtdf.slc.tmatch.idx >= 0) & (selected_evtdf.mccoh)
bench = "nu coherent"
bench_df = selected_evtdf[is_coh]
print_df(bench_df,bench)

is_cosmic = (selected_evtdf.slc.tmatch.idx < 0)
bench = "Cosmic"
bench_df = selected_evtdf[is_cosmic]
print_df(bench_df,bench)

