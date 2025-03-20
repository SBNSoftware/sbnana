# How to run this script:

# $ python unc_run_evtSel.py [selection_keys] CL pc_uncertainty_sig pc_uncertainty_bg
#selection_key is just the letter representing the key in evtSel_dict, defined in unc_cuts.py
#CL is confidence level passed as a decimal (default is 0.9)
#pc_uncertainty_sig is total signal uncertainty passed as a decimal (default is 0.5)
#pc_uncertainty_bg is total background uncertainty passed as a decimal (default is 0.5)


# ------------------------------------------------------------------------

import uproot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
import h5py
import math
import subprocess
import sys
import pyhf
import re
import io

from util import *

import var
import cut
import data
import hist

import importlib

from pyanalib import panda_helpers
from unc_funcs import *
from unc_samples import *
from unc_MC_overhead import *
from unc_other_limits import *
from unc_cuts import *

# ------------------------------------------------------------------------
# HPS Specifics

#mass_to_hpsSampleIndex = dict(zip(higgs_masses, range(len(higgs_masses))))
mass_to_hpsSampleIndex = dict(zip(np.array(higgs_masses)/1000., range(len(higgs_masses)))) # evtdf or mchdf
def expected_hps_events(selected_evtdf, mass, new_theta): # HPS
    categories = make_categories(selected_evtdf, detailed_bsm=True)
    i = mass_to_hpsSampleIndex[mass]
    df = selected_evtdf[categories[i]]
    old_th = float(higgs_mcdfs[i].iloc[[0]].C1)
    rescale_new_mixing = []
    for idx in df.index:
        row = higgs_mcdfs[i].loc[(idx[1],idx[2], 0)]
        factor = float(reweight_mixing(new_theta, row.start, row.enter, row.exit, row.decay_length, old_th))
        rescale_new_mixing.append(factor)
    x = sum(np.array(df.scale)*np.array(rescale_new_mixing))
    return(x)

# ------------------------------------------------------------------------
# ALP Specifics

mass_to_alpSampleIndex = dict(zip(
    np.array(alp_nosup_masses)/1000.,
    np.array(range(len(alp_nosup_masses)))+np.array([len(higgs_masses)]*len(alp_nosup_masses)) 
))
mass_to_alp_mcdf_SampleIndex = dict(zip(np.array(alp_nosup_masses)/1000., range(len(alp_nosup_masses)))) # mchdf

# get the X,Y pairs for running coupling for cmu with fa:
running_cmu_codominance_file = '/exp/icarus/data/users/jdyer/muon_selection/tabulated_running_alp_cmu_vals_fromJosh/cl_c3_1_c12_1.csv' 
Y_alp = [] # fa values
fa_strings = []
param_pair_strings = []
fa_and_running_cmu_pairs_codominance = []
with open(running_cmu_codominance_file, 'r') as file:
    for line in file: 
        fa, cl = line.strip().split(',')
        fa_and_running_cmu_pairs_codominance.append( (float(fa), float(cl)) )
        Y_alp.append(float(fa))
        fa_strings.append('fa='+str(fa))
        param_pair_strings.append('1/f_a='+str(1./float(fa))+', c_mu='+str(cl))
X_alp_fromP0 = np.array(alp_nosup_masses).astype(float)/1000. # alp masses to scan over, in GeV
X_alp_fromK = np.array(higgs_masses).astype(float)/1000. # alp masses to scan over, in GeV
        
def expected_alp_events(selected_evtdf, mass, new_fa, new_cl, print_stuff=False, do_P0prod=True, do_Kprod=True): # input mass in GeV
    categories = make_categories(selected_evtdf, detailed_bsm=True)

    # REWEIGHT FROM MESON-MIXING-PRODUCED ALPs: 
    
    if do_P0prod:
        if mass in mass_to_alpSampleIndex.keys():
            i_alp = mass_to_alpSampleIndex[mass]
            i_alp_mcdf = mass_to_alp_mcdf_SampleIndex[mass]
            df_alp = selected_evtdf[categories[i_alp]]
            old_fa = float(alp_nosup_mcdfs[i_alp_mcdf].iloc[[0]].C1)
            old_cmu = float(alp_nosup_mcdfs[i_alp_mcdf].iloc[[0]].C2)        
            old_alp_f = float(alp_nosup_mcdfs[i_alp_mcdf].iloc[0].allowed_decay_fraction)    
            reweigts_psuedomixed_alps = []
            for idx in df_alp.index:
                row = alp_nosup_mcdfs[i_alp_mcdf].loc[(idx[1],idx[2], 0)]
                factor_rwgt_producedAlps = float(
                    reweight_alps(old_fa, new_fa, old_cmu, new_cl, 
                                  row.start, row.enter, row.exit, 
                                  row.decay_length, old_alp_f)
                )
                reweigts_psuedomixed_alps.append(factor_rwgt_producedAlps)
            alps_from_mixing = sum(np.array(df_alp.scale)*np.array(reweigts_psuedomixed_alps))
        else: 
            print('We cannot perform a reweighting from ALPs -> ALPs for that mass, because we did not generate meson-mixing-produced ALPs for that mass.')
            alps_from_mixing = -1.
    else: alps_from_mixing = -2.
        
    # REWEIGHT FROM KAON-DECAY-PRODUCED ALPs:
    
    if do_Kprod: 
        if mass < klong_mass-pizero_mass:
            if mass in mass_to_hpsSampleIndex.keys():
                i_hps = mass_to_hpsSampleIndex[mass]
                df_hps = selected_evtdf[categories[i_hps]]
                df_hps = df_hps[np.abs(df_hps.slc.truth.parent_pdg == 321)] ## Mask out the ones with K0_L parent since that is CP-violating-suppressed for ALPs. # This line implemented Feb. 3, 2025.

                old_theta = float(higgs_mcdfs[i_hps].iloc[[0]].C1)
                hps_f = float( HPS_MuonPartialWidth( float(higgs_mcdfs[i_hps].iloc[[0]].M), float(higgs_mcdfs[i_hps].iloc[[0]].C1))/
                              HPS_total_decay_Width(float(higgs_mcdfs[i_hps].iloc[[0]].M), float(higgs_mcdfs[i_hps].iloc[[0]].C1)), 
                             )
                if mass in mass_to_alpSampleIndex.keys(): 
                    i_alp_mcdf = mass_to_alp_mcdf_SampleIndex[mass]
                    alp_f = float(alp_nosup_mcdfs[i_alp_mcdf].iloc[0].allowed_decay_fraction)
                else:
                    total_a_to_uu = ALP_total_decay_Width(mass, new_fa, new_cl, 1,1,1)
                    assert total_a_to_uu != -1
                    alp_f =ALP_uu_width(mass, new_fa, new_cl)/ALP_total_decay_Width(mass, new_fa, new_cl, 1,1,1)
                reweigts_kdecay_alps = []
                for idx in df_hps.index:
                    row = higgs_mcdfs[i_hps].loc[(idx[1],idx[2], 0)]
                    hps_mean_dist = row.decay_length
                    factor_rwgt_hps = float(
                        reweight_hps_to_alps(
                            mass, old_theta, hps_mean_dist, hps_f, 
                            new_fa, new_cl, alp_f,
                            row.start, row.enter, row.exit, int(df_hps.loc[idx].slc.truth.parent_pdg)
                        )
                    )
                    reweigts_kdecay_alps.append(factor_rwgt_hps)
                alps_from_kdecay = sum(np.array(df_hps.scale)*np.array(reweigts_kdecay_alps))
            else:
                print('input mass: %a' % mass, klong_mass-pizero_mass)
                print('Production via kaon decay IS possible, but we did not generate HPSs at the desired alp mass, so we cannot perform the needed rescaling.' )
                alps_from_kdecay = -1.
        else:
            alps_from_kdecay = 0.
    else: alps_from_kdecay = -2.
                  
    # Combine both Production Modes and return results:
    
    if (alps_from_mixing == -1.) | (alps_from_kdecay == -1.):
        alps_from_both = -1.
    elif (alps_from_mixing == -2.) | (alps_from_kdecay == -2.):
        alps_from_both = -2.
    else:
        alps_from_both = alps_from_mixing + alps_from_kdecay
                  
    return(alps_from_mixing, alps_from_kdecay, alps_from_both)

# ------------------------------------------------------------------------
# pyhf
def return_modelspec(mc_signal, pc_uncertainty_sig, total_selected_bg, pc_uncertainty_bg):
    modelspec = {
    "channels": [
        {
        "name": "singlechannel",
        "samples": [
            {
            "name": "signal",
            "data": [ mc_signal ],
            "modifiers": [
                {"name": "mu", "type": "normfactor", "data": None },
                {"name": "uncorr_siguncrt", "type": "shapesys", "data": [mc_signal*pc_uncertainty_sig]}
                ]
            },
            {
                "name": "background",
                "data": [ total_selected_bg ],       # total_selected_bg
                "modifiers": [
                { "name": "uncorr_bkguncrt", "type": "shapesys", "data": [total_selected_bg*pc_uncertainty_bg] }
                                     # if equal to data, then 100% uncertainty.

                ]
          }
        ]
      }
    ]
    }
    return modelspec

# ------------------------------------------------------------------------
# Main function:
def pushEventSel(evtdf, selection_key, CL=0.9, pc_uncertainty_sig=0.5, pc_uncertainty_bg=0.5):#, obs_data=2):
    # selection_key should be a string - it's the key to evtSel_dict defined in unc_cuts.py

    selection_name = selection_key
    cut_list = evtSel_dict[selection_key][0]
    thresholds = evtSel_dict[selection_key][1]
    selection_has_been_run = False
    output_path = '/exp/icarus/data/users/jdyer/muons_selections_study_2501/Selection_'+selection_name+'/'
    try:
        result = subprocess.run(['ls', output_path], capture_output=True, text=True, check=True)
        selection_has_been_run = True
        #print(result.stdout)
        selected_evtdf = pd.read_pickle(output_path + 'selected_evtdf')
        res_pot = pd.read_pickle(output_path + 'res_pot')
        res_mc = pd.read_pickle(output_path + 'res_mc')
        res_pc = pd.read_pickle(output_path + 'res_percent')

    except subprocess.CalledProcessError as e:
        subprocess.run(['mkdir', output_path])
        #print(f"Error: {e}", f"Return code: {e.returncode}", f"Output: {e.stdout}", f"Error output: {e.stderr}", sep='\n')

        # Impose that one or both tracks are uncontained:
        when_uncontained = ~TrkInFV(evtdf.trunk.trk.end) | ~TrkInFV(evtdf.branch.trk.end)
        evtdf = evtdf[when_uncontained]

        # Apply the event selection
        cut_results = apply_cuts(evtdf, cut_list, thresholds=thresholds, detailed_nu='none')
        mask = cut_results[-1]
        res_mc = cut_results[0]
        res_pot = cut_results[1]
        res_pc = cut_results[2]
        selected_evtdf = add_hdr_info(evtdf[mask], hdrs)

        selected_evtdf.to_pickle(output_path + 'selected_evtdf')
        res_mc.to_pickle(output_path + 'res_mc')
        res_pot.to_pickle(output_path + 'res_pot')
        res_pc.to_pickle(output_path + 'res_percent')
          
# Model-Independent Sensitivity Stuff:
    # depends on CL, and uncertainties, but not models.
    # Does need total_selected_bg, which is event selection dependent.
    # Note: in real analysis, the signal uncertainties for HPS and ALP will probs differ.
    
    cats = make_categories(selected_evtdf)
    total_selected_bg = np.sum(selected_evtdf[cats[1]].scale) + np.sum(selected_evtdf[cats[2]].scale) # selected nus and cosmics
    observed_data_events = math.ceil(total_selected_bg) # pretend for now
    mc_signal = 11 # should be in ballpark of number of signal events you expect to exclude with the desired exclusion limit.
    confidence = 1-CL # Desired Confidence Level: The confidence level is 1 minus this value.
    limitname = "pyhf_ExpLim_%a CL_ %a uncSig_ %a uncBg" % (int(CL*100), int(pc_uncertainty_sig*100), int(pc_uncertainty_bg*100))
    limitname = limitname.replace(' ','')
    try:
        exp_limits = np.load(output_path+limitname+'.npy')
    except:
        model = pyhf.Model( return_modelspec(mc_signal, pc_uncertainty_sig, total_selected_bg, pc_uncertainty_bg) )
        poi_values = np.linspace(0.1, 10.0, 200)
        obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(
            [observed_data_events] + model.config.auxdata,
            model,
            poi_values,
            level=confidence,
            return_results=True
        )
        np.save(output_path+limitname, exp_limits)

# HPS SENSITIVITY
    hps_plotname = "HPS_"+selection_key+"_ %a CL_ %a uncSig_ %a uncBg.png" % (int(CL*100), int(pc_uncertainty_sig*100), int(pc_uncertainty_bg*100))
    hps_plotname = hps_plotname.replace(' ','')
    try:
        result = subprocess.run(["ls", output_path+hps_plotname], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        # Find new expected number of events as function of model parameters:
        X = np.array(higgs_masses).astype(float)
        Y = np.logspace(-6,-1,40)
        try: # if selection_has_been_run:
            Z = np.load(output_path+'hps_Z.npy')
        except: # else:
            Z = np.reshape([expected_hps_events(selected_evtdf, m/1000., theta) for theta in Y  for m in X ],(40,8))
            np.save(output_path+'hps_Z', Z)
        
        # Make plot
        plt.plot(hps_na62_x, hps_na62_y, label='NA62', color='purple')
        plt.plot(hps_uB_x, hps_uB_y, label='uBooNE', color='#FFDB58')
        plt.plot(hps_E949_x, hps_E949_y, label='E949', color='red')
        plt.plot(hps_LHCb_x, hps_LHCb_y, label='LHCb', color='C1')
        plt.plot(hps_Gray_x, hps_Gray_y, label='ICARUS cont. $\mu\mu$ search', color='C0')
        mycolor = 'black'
        plt.contourf(X,Y*Y,np.log(Z),levels=[np.log(exp_limits[0]*mc_signal),np.log(exp_limits[4]*mc_signal)],colors=mycolor, alpha=0.2)
        plt.contourf(X,Y*Y,np.log(Z),levels=[np.log(exp_limits[1]*mc_signal),np.log(exp_limits[3]*mc_signal)],colors=mycolor, alpha=0.4)
        plt.contour(X,Y*Y,np.log(Z),levels=[np.log(exp_limits[2]*mc_signal)],colors=mycolor)
        plt.plot([500, 600],[1,2], label='ICARUS w/ uncont. $\mu\mu$ \nexpectation $\pm 1,2\sigma$', color=mycolor)
        
        # General formatting
        plt.title('Selection '+selection_name+'\n%a%% CL \nassuming $\\delta_{sig}$=%a%%, $\\delta_{bg}$=%a%%' % (CL*100, pc_uncertainty_sig*100, pc_uncertainty_bg*100) )
        plt.yscale('log')
        plt.xlabel('$m_S$ (MeV)')
        plt.ylabel('$\\theta^2$')#_S$')
        plt.legend()
        plt.xlim((220, 340))
        plt.ylim((1e-10,2e-6))
        plt.savefig(output_path + hps_plotname, format='png', bbox_inches='tight')
        plt.close()
          
# ALP SENSITIVITY (assuming running coupling between cmu and fa.)

    alp_plotname = "ALP_"+selection_key+"_ %a CL_ %a uncSig_ %a uncBg.png" % (int(CL*100), int(pc_uncertainty_sig*100), int(pc_uncertainty_bg*100))
    alp_plotname = alp_plotname.replace(' ','')
    try:
        result = subprocess.run(["ls", output_path+alp_plotname], capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        # Find new expected number of events as function of model parameters:
        
        # From P0:
        try:
            Z_alp_fromP0 = np.load(output_path+'alp_Z_fromP0')
        except: # else:
            Z_alp_fromP0 = [expected_alp_events(selected_evtdf, m, pair[0], pair[1], do_Kprod=False) for pair in fa_and_running_cmu_pairs_codominance for m in X_alp_fromP0]
            Z_alp_fromP0 = np.array(Z_alp_fromP0)
            np.save(output_path+'alp_Z_fromP0', Z_alp_fromP0)
        Z_inds_fromP0 = np.arange(len(Z_alp_fromP0))
        Z_outputs_per_mass_fromP0 = []
        Z_mixProd_fromP0 = []
        for nm in range(len(X_alp_fromP0)):
            Z_of_this_mass = Z_alp_fromP0[ np.where(Z_inds_fromP0%len(X_alp_fromP0) == nm) ]
            Z_outputs_per_mass_fromP0.append( Z_of_this_mass )
            Z_mixProd_fromP0.append(Z_of_this_mass[:,0])
        Z_mixProd_fromP0 = np.transpose( np.array(Z_mixProd_fromP0) )
   

        # From K: 
        try:
            Z_alp_fromK = np.load(output_path+'alp_Z_fromK')
        except:
            Z_alp_fromK = [expected_alp_events(selected_evtdf, m, pair[0], pair[1], do_P0prod=False) for pair in fa_and_running_cmu_pairs_codominance for m in X_alp_fromK]
            Z_alp_fromK = np.array(Z_alp_fromK)
            np.save(output_path+'alp_Z_fromK', Z_alp_fromK)
            
        Z_inds_fromK = np.arange(len(Z_alp_fromK))
        Z_outputs_per_mass_fromK = []
        Z_KdecayProd_fromK = []
        for nm in range(len(X_alp_fromK)):
            Z_of_this_mass = Z_alp_fromK[ np.where(Z_inds_fromK%len(X_alp_fromK) == nm) ]
            Z_outputs_per_mass_fromK.append( Z_of_this_mass )
            Z_KdecayProd_fromK.append(Z_of_this_mass[:,1])
        Z_KdecayProd_fromK = np.transpose( np.array(Z_KdecayProd_fromK) )

        # Make plot:
        
        fig = plt.figure()
        ax = plt.subplot(111)
        
        # Vertical lines:
        myalpha = 0.5
        mylinestyle = ':'
        vcolor = 'brown'
        ## plt.axvline([klong_mass-pizero_mass], linestyle=mylinestyle, color=vcolor, linewidth=2, alpha=myalpha) # where kaon-decay production mode becomes impossible. # Nix -- this mode is CP suppressed.
        plt.axvline([kplus_mass-piplus_mass], linestyle=mylinestyle, color=vcolor, linewidth=2, alpha=myalpha) # where kaon-decay production mode becomes impossible.
        plt.axvline([0.548], linestyle=mylinestyle, color=vcolor, linewidth=2, alpha=myalpha) # eta mass
        plt.axvline([3*pizero_mass], linestyle=mylinestyle, color=vcolor, linewidth=2, alpha=myalpha) # axion mass where hadronic decays start to matter

        # Other limits:
        plt.plot(alp_NA62_x, alp_NA62_y, color='purple', label='NA62') # 
        plt.plot(alp_uB_x, alp_uB_y, color='#FFDB58', label='uBooNE') # 
        plt.plot(alp_charm_uu_x, alp_charm_uu_y, color='green', label='CHARM uu, from K')
        plt.plot(alp_charm_gg_x, alp_charm_gg_y, color='cyan', label='CHARM gg, from K')
        #plt.plot(charmP0_uu_x, charmP0_uu_y, color='orange', label='CHARM uu, from P0')
        plt.plot(charmP0_uu_x_1, charmP0_uu_y_1, color='orange', label='CHARM uu, from P0')
        plt.plot(charmP0_uu_x_2, charmP0_uu_y_2, color='orange')
        #plt.plot(charmP0_gg_x, charmP0_gg_y, color='red', label='CHARM gg, from P0')
        plt.plot(charmP0_gg_x_1, charmP0_gg_y_1, color='red', label='CHARM gg, from P0')
        plt.plot(charmP0_gg_x_2, charmP0_gg_y_2, color='red')
        plt.plot(alp_Gray_x/1000., alp_Gray_y, label='ICARUS cont. $\mu\mu$ search', color='C0')
        
        # My limits:
        mycolor = 'black'
        
        # from mixing
        plt.contourf(X_alp_fromP0, 1./np.array(Y_alp), np.log(Z_mixProd_fromP0), levels=[np.log(exp_limits[0]*mc_signal),np.log(exp_limits[4]*mc_signal)],colors=mycolor, alpha=0.2)
        plt.contourf(X_alp_fromP0, 1./np.array(Y_alp), np.log(Z_mixProd_fromP0), levels=[np.log(exp_limits[1]*mc_signal),np.log(exp_limits[3]*mc_signal)],colors=mycolor, alpha=0.4)
        plt.contour(X_alp_fromP0, 1./np.array(Y_alp), np.log(Z_mixProd_fromP0), levels=[np.log(exp_limits[2]*mc_signal)],colors=mycolor, linestyles='dotted')
        plt.plot([X_alp_fromP0[0], X_alp_fromP0[1]],[1,2], label="This analysis, from mixing", color=mycolor, linestyle=':', alpha=1)

        # from K-decay
        plt.contourf(X_alp_fromK, 1./np.array(Y_alp), np.log(Z_KdecayProd_fromK), levels=[np.log(exp_limits[0]*mc_signal),np.log(exp_limits[4]*mc_signal)],colors=mycolor, alpha=0.2)
        plt.contourf(X_alp_fromK, 1./np.array(Y_alp), np.log(Z_KdecayProd_fromK), levels=[np.log(exp_limits[1]*mc_signal),np.log(exp_limits[3]*mc_signal)],colors=mycolor, alpha=0.4)
        plt.contour(X_alp_fromK, 1./np.array(Y_alp), np.log(Z_KdecayProd_fromK), levels=[np.log(exp_limits[2]*mc_signal)],colors=mycolor, linestyles='dashed')
        plt.plot([X_alp_fromK[0], X_alp_fromK[1]],[1,2], label="This analysis, from K-decay", color='black', linestyle='--', alpha=1)
        
        # general formatting
        plt.title('Selection '+selection_name+'\n%a%% CL \nassuming $\\delta_{sig}$=%a%%, $\\delta_{bg}$=%a%%' % (CL*100, pc_uncertainty_sig*100, pc_uncertainty_bg*100) )
        plt.yscale('log')
        plt.xlabel('$m_S$ (GeV)')
        plt.ylabel('$1/f_a$')
        plt.xlim((0.22, 0.45))
        plt.ylim((1e-6, max(1./np.array(Y_alp))))
        #plt.legend(loc='lower right')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.2), ncol=3, fontsize='small')
        plt.savefig(output_path + alp_plotname, format='png', bbox_inches='tight')
        plt.close()

# end of main function.
        
selection_keys = sys.argv[1]
output = io.StringIO() # Create a string buffer
print(re.findall(r"\w+",selection_keys), file=output) 
    # Redirect the output of print to the buffer
    # Note: r"\w+" matches one or more consecutive alphanumeric characters (letters, digits, and underscores).
selection_keys = output.getvalue() # Get the string from the buffer
output.close() # Close the buffer
selection_keys = eval(selection_keys)

#str_to_eval = 're.findall(r"\w+","'+selection_keys+'")'
#selection_keys = eval('re.findall(r"\w+", '+selection_keys+')')
if __name__ == "__main__":
    for key in selection_keys:
        if len(sys.argv) > 2:
            CL = float(sys.argv[2])
            pc_uncertainty_sig = float(sys.argv[3])
            pc_uncertainty_bg = float(sys.argv[4])
            pushEventSel(evtdf, key, CL=CL, pc_uncertainty_sig=pc_uncertainty_sig, pc_uncertainty_bg=pc_uncertainty_bg)
        else:
            pushEventSel(evtdf, key)

            
            
            
            