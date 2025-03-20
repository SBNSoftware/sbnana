# Useful for: 
#    unc_VeryFarSB_McData_w_cohlike.ipynb

# Here, I use the custom function 'data.mc_dataset(mccohfile, "evt", mcnukey="mcnuweight").df' to get the evt dataframes from the files listed in unc_samples.py. 

# The function takes a while, because it uses information in [] to compute weights needed to run the systematics, and then inserts those weights into the returned dataframe. 

# The returned dataframes are also very large, and it is not practical to handle so much information for the complete dataset.

# So, I have the following workflow: ONE TIME, I load the dataframes from their original files with the mc_dataset function, and then I restrict that dataset to the Far Sideband to shrink it considerable and make it handle-able in the notebook. I save the resulting dataframes at:

# []

# All subsequent times I run the notebook, I load in those systematics/FarSB dataframes instead, saving time. 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SOME GENERAL OVERHEAD STUFF:

## test
#f = h5py.File(mccohfile, 'r')
#print(f.keys())

import pandas as pd
import pickle
import math
import var
import data
from unc_samples import *
from unc_funcs import *

NuMI_angle_thresh = 15 # deg # Very Far Sideband I presented in 2023: > 0.15 rad # Use for selection: < 0.05 rad
open_angle_thresh = [70,115]# Use for selection: < 0.35 rad 

def sb_numi_angle_mask(df, thresh = NuMI_angle_thresh): # do calculation in radians. User uses degrees.
    return df.Snumi_angle_mcs > thresh*math.pi/180., 'S_NuMI_angle > '+str(thresh)+' deg'#'\u00B0'

def sb_open_angle_mask(df, thresh = open_angle_thresh):
    return ( ( (np.arccos(dotdf(df.trunk.trk.dir, df.branch.trk.dir)) > thresh[0]*math.pi/180.) & 
            (np.arccos(dotdf(df.trunk.trk.dir, df.branch.trk.dir)) < thresh[1]*math.pi/180.) ), 
            str(thresh[0])+' deg < opening angle < '+str(thresh[1])+' deg')

def veryFarSB(df):
    df = df.copy()
    print("At start: ", df.shape)
    sideband_mask = sb_numi_angle_mask(df)[0] & sb_open_angle_mask(df)[0]
    df = df[sideband_mask]
    print("After restriction to sideband: ", df.shape)
    return df

def resave_datasets(df_names, df_list, output):
    with pd.HDFStore(output) as hdf:
        for k,df in zip(df_names, df_list):
            try:
                hdf.put(key=k, value=df, format="fixed")
            except Exception as e:
                print("Table %s failed to save, skipping. Exception: %s" % (k, str(e)))
    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HERE IS THE STUFF TO RUN ONCE:
# Note 24/09/17: I have successfully tested this for mccoh sample (and it looked okay) since cleaning up/moving to this overhead file, but not for the others due to storage issues. The samples I made previously have been moved to: /exp/icarus/data/users/jdyer/muon_datamc/samples_240909/

#savedir = "/exp/icarus/data/users/jdyer/muon_datamc/"
savedir_sb = "/exp/icarus/data/users/jdyer/muon_datamc/lightened_wConcreteFix/"
savedir_all_mc = "/exp/icarus/data/users/jdyer/muon_datamc/allMC_ie_no_restriction_to_sideband_wConcreteFix/"


# ~~~~~~
## MC Incoh:
## takes about 23 minutes to run.

#mc_incoh_evt = data.mc_dataset(nu_file, "evt", mcnukey="mcnuwgt", mccut=~is_coh_like, mccut_any=False)
#mc_incoh_evt_df = mc_incoh_evt.df
#print('mc_incoh livetime: ', mc_incoh_evt.livetime)
#print('mc_incoh POT: ', mc_incoh_evt.POT)
#print(mc_incoh_evt_df.shape)
## Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events):
#mc_incoh_evt_df = mc_incoh_evt_df[satisfies_new_FV(mc_incoh_evt_df)] 
## Require exiting:
#when_uncontained = ~TrkInFV(mc_incoh_evt_df.trunk.trk.end) | ~TrkInFV(mc_incoh_evt_df.branch.trk.end)
#mc_incoh_evt_df = mc_incoh_evt_df[when_uncontained]
#mc_incoh_evt_df_all = add_calculated_evtdf_cols(mc_incoh_evt_df)#, newcol_name=None, newcol_val=None)
#mc_incoh_evt_df_sb = veryFarSB(mc_incoh_evt_df_all)

##lightened_mc_incoh_evt_df.to_pickle(savedir + 'new_mc_incoh_evt_df') # intermediate pickle save.
##all_mc_incoh_evt_df.to_pickle(savedir +'allMC_ie_no_restriction_to_sideband/new_mc_incoh_evt_df')
### If I've already made the dataframe and pickled it, then read it in here:
##lightened_mc_incoh_evt_df = pd.read_pickle(savedir + 'new_mc_incoh_evt_df')

##mcdf = pd.read_hdf(nu_file, key='mcnuwgt')
#mc_incoh_hdr_df = pd.read_hdf(nu_file, "hdr")
#mc_incoh_mcnu_df = pd.read_hdf(nu_file, "mcnu")

#df_names = ['evt', 'hdr', 'mcnu']
#resave_datasets(df_names, 
#                [mc_incoh_evt_df_all, mc_incoh_hdr_df, mc_incoh_mcnu_df],
#                savedir_all_mc+'mc_incoh_all'
#               )
#resave_datasets(df_names, 
#                [mc_incoh_evt_df_sb, mc_incoh_hdr_df, mc_incoh_mcnu_df],
#                savedir_all_mc+'mc_incoh_sb'
#               )

#del mc_incoh_evt
#del mc_incoh_evt_df
#del mc_incoh_evt_df_all
#del mc_incoh_evt_df_sb

# ~~~~~~
## MC Coh:
## takes about 7 minutes to run.

## Notes: 
#Keys I do: evt, hdr, mch, mcnu

#This sample contains coh-LIKE events. Really, it is any event with a muon or a pion in the final state. This means that the interaction isn't neccessarily a coherent interaction, but it just looks like one. 

#Note: It is possible that there is pileup in the neutrino interactions. If it happens that two neutrino interactions overlap in a single reconstructed slice, then even if the "bonus" one is not coh-like, it would be included in this sample as well. Having two neutrino interactions overlap like this is extremely unlikely, so this probably never happens.

#mccohfile_evt = data.mc_dataset(cohlike_nu_file, "evt", mcnukey="mcnuwgt")
#mccohfile_evt_df = mccohfile_evt.df
#print('mccoh livetime: ', mccohfile_evt.livetime)
#print('mccoh POT: ', mccohfile_evt.POT)
#print(mccohfile_evt_df.shape)
## Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events):
#mccohfile_evt_df = mccohfile_evt_df[satisfies_new_FV(mccohfile_evt_df)] 
## Require exiting:
#when_uncontained = ~TrkInFV(mccohfile_evt_df.trunk.trk.end) | ~TrkInFV(mccohfile_evt_df.branch.trk.end)
#mccohfile_evt_df = mccohfile_evt_df[when_uncontained]
#mccohfile_evt_df_all = add_calculated_evtdf_cols(mccohfile_evt_df)#, newcol_name=None, newcol_val=None)
#mccohfile_evt_df_sb = veryFarSB(mccohfile_evt_df_all)

#mccohfile_hdr_df = pd.read_hdf(cohlike_nu_file, "hdr")
#mccohfile_mch_df = pd.read_hdf(cohlike_nu_file, "mch")
#mccohfile_mcnu_df = pd.read_hdf(cohlike_nu_file, "mcnu")

## resave dataset with lightened evt dataframe:
#df_names = ['evt', 'hdr', 'mch', 'mcnu']
#resave_datasets(df_names, 
#                [mccohfile_evt_df_all, mccohfile_hdr_df, mccohfile_mch_df, mccohfile_mcnu_df], 
#                savedir_all_mc+'mc_coh_all'
#               )
#resave_datasets(df_names, 
#                [mccohfile_evt_df_sb, mccohfile_hdr_df, mccohfile_mch_df, mccohfile_mcnu_df], 
#                savedir_sb+'mc_coh_sb'
#               )
#del mccohfile_evt
#del mccohfile_evt_df
#del mccohfile_evt_df_all
#del mccohfile_evt_df_sb

# ~~~~~~
## MC Coh Detector Variation Samples:

#for i, file in enumerate(cohlike_nu_detVar_files):
#    print(file, syst_labels[i], sep=' ')
#    if ( i==0 | i ==1 | i==6): continue
#    #f = h5py.File(systFile, 'r')
#    #print(f.keys(), '\n')
    
#    evt = data.mc_dataset(file, "evt", mcnukey="mcnuwgt")
#    print('\n')
#    print('systFile: ', file)
#    print('livetime: ', evt.livetime)
#    print('POT: ', evt.POT)
#    evt_df = evt.df

#    # Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events):
#    evt_df = evt_df[satisfies_new_FV(evt_df)] 
#    # Require exiting:
#    when_uncontained = ~TrkInFV(evt_df.trunk.trk.end) | ~TrkInFV(evt_df.branch.trk.end)
#    evt_df = evt_df[when_uncontained]

#    evt_df = add_calculated_evtdf_cols(evt_df)
#    evt_df_sb = veryFarSB(evt_df)
#    #evt_df_sb.to_pickle(savedir + 'syst_mccoh_evt_df_'+syst_labels[i].replace(' ', '_'))
    
#    # resave:
#    output_sb = savedir_sb+detVar_labels[i].replace(' ', '_')
#    output_all_mc = savedir_all_mc+detVar_labels[i].replace(' ', '_')
#    for output in [output_sb, output_all_mc]:
#        try: 
#            subprocess.run(['rm', output])
#        except:
#            print("Did not delete a pre-existing file.")
#        with pd.HDFStore(output) as hdf:
#            # save the evt key to the dataset:
#            if output == output_sb:
#                evtdf = evt_df_sb
#            else: 
#                evtdf = evt_df
#            try:
#                hdf.put(key='evt', value=evtdf, format="fixed")
#            except Exception as e:
#                print("Table %s failed to save, skipping. Exception: %s" % ('evt', str(e)))
#            # save the hdr key to the dataset:    
#            try:
#                hdf.put(key='hdr', value=pd.read_hdf(file, "hdr"), format="fixed")
#                hdf.put(key='mcnu', value=pd.read_hdf(file, "mcnu"), format="fixed")
#            except Exception as e:
#                print("Table %s failed to save, skipping. Exception: %s" % ('hdr', str(e)))
        
# ~~~~~~
## Onbeam Data:  

##f = h5py.File(onbeamfile_Run1, 'r')
##print(f.keys(), '\n')

#onbeamR1_evt = data.onbeam_dataset(onbeamfile_Run1, "evt")
#print('onbeamR1 POT: ', onbeamR1_evt.POT)
#print('onbeamR1 livetime: ', onbeamR1_evt.livetime)
#onbeamR1_evt_df = onbeamR1_evt.df
#onbeamR1_evt_df_sb = veryFarSB(onbeamR1_evt_df)
## Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events):
#onbeamR1_evt_df_sb = onbeamR1_evt_df_sb[satisfies_new_FV(onbeamR1_evt_df_sb)] 
## Require exiting:
#when_uncontained = ~TrkInFV(onbeamR1_evt_df_sb.trunk.trk.end) | ~TrkInFV(onbeamR1_evt_df_sb.branch.trk.end)
#onbeamR1_evt_df_sb = onbeamR1_evt_df_sb[when_uncontained]
#onbeamR1_evt_df_sb = add_calculated_evtdf_cols(onbeamR1_evt_df_sb, 
#                                               newcol_name=["Run1", "Run2"], newcol_val=[True, False])

#onbeamR2_evt = data.onbeam_dataset(onbeamfile_Run2, "evt")
#print('onbeamR2 POT: ', onbeamR2_evt.POT)
#print('onbeamR2 livetime: ', onbeamR2_evt.livetime)
#onbeamR2_evt_df = onbeamR2_evt.df
#onbeamR2_evt_df_sb = veryFarSB(onbeamR2_evt_df)
## Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events):
#onbeamR2_evt_df_sb = onbeamR2_evt_df_sb[satisfies_new_FV(onbeamR2_evt_df_sb)] 
## Require exiting:
#when_uncontained = ~TrkInFV(onbeamR2_evt_df_sb.trunk.trk.end) | ~TrkInFV(onbeamR2_evt_df_sb.branch.trk.end)
#onbeamR2_evt_df_sb = onbeamR2_evt_df_sb[when_uncontained]

#onbeamR2_evt_df_sb = add_calculated_evtdf_cols(onbeamR2_evt_df_sb, 
#                                               newcol_name=["Run1", "Run2"], newcol_val=[False, True]))

## resave:
#for run_string in ['offbeamR1', 'offbeamR2']:
#    output = savedir_sb+run_string
#    try: 
#        subprocess.run(['rm', output])
#    except:
#        print("Did not delete a pre-existing file.")
#    with pd.HDFStore(output) as hdf:
#        if run_string == 'offbeamR1':
#            evtdf = 
#            hdr = pd.read_hdf(onbeamfile_Run1, "hdr")
#        if run_string == 'offbeamR2': 
#            evtdf = onbeamR2_evt_df_sb
#            hdr = pd.read_hdf(onbeamfile_Run2, "hdr")
#        try:
#            hdf.put(key='evt', value=evtdf, format="fixed")
#            hdf.put(key='hdr', value=hdr, format="fixed")
#        except Exception as e:
#            print("Table %s failed to save, skipping. Exception: %s" % ('hdr', str(e)))

#del onbeamR1_evt
#del onbeamR1_evt_df
#del onbeamR1_evt_df_sb

#del onbeamR2_evt
#del onbeamR2_evt_df
#del onbeamR2_evt_df_sb

# ~~~~~~
## Offbeam Data: 

#offbeamR1_evt = data.offbeam_dataset(offbeamfile_Run1, "evt")
#print('offbeamR1 livetime: ', offbeamR1_evt.livetime)
#offbeamR1_evt_df = offbeamR1_evt.df
#offbeamR1_evt_df_sb = veryFarSB(offbeamR1_evt_df)
## Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events):
#offbeamR1_evt_df_sb = offbeamR1_evt_df_sb[satisfies_new_FV(offbeamR1_evt_df_sb)] 
## Require exiting:
#when_uncontained = ~TrkInFV(offbeamR1_evt_df_sb.trunk.trk.end) | ~TrkInFV(offbeamR1_evt_df_sb.branch.trk.end)
#offbeamR1_evt_df_sb = offbeamR1_evt_df_sb[when_uncontained]
#offbeamR1_evt_df_sb = add_calculated_evtdf_cols(offbeamR1_evt_df_sb, 
#                                               newcol_name=["Run1", "Run2"], newcol_val=[True, False])

#offbeamR2_evt = data.offbeam_dataset(offbeamfile_Run2, "evt")
#print('offbeamR2 livetime: ', offbeamR2_evt.livetime)
#offbeamR2_evt_df = offbeamR2_evt.df
#offbeamR2_evt_df_sb = veryFarSB(offbeamR2_evt_df)
## Make sure FV definition is up to date (so that we don't accept any too-short non-fiducial events):
#offbeamR2_evt_df_sb = offbeamR2_evt_df_sb[satisfies_new_FV(offbeamR2_evt_df_sb)] 
## Require exiting:
#when_uncontained = ~TrkInFV(offbeamR2_evt_df_sb.trunk.trk.end) | ~TrkInFV(offbeamR2_evt_df_sb.branch.trk.end)
#offbeamR2_evt_df_sb = offbeamR2_evt_df_sb[when_uncontained]
#offbeamR2_evt_df_sb = add_calculated_evtdf_cols(offbeamR2_evt_df_sb, 
#                                               newcol_name=["Run1", "Run2"], newcol_val=[False, True])

## resave:
#for run_string in ['onbeamR1', 'onbeamR2']:
#    output = savedir_sb+run_string
#    try: 
#        subprocess.run(['rm', output])
#    except:
#        print("Did not delete a pre-existing file.")
#    with pd.HDFStore(output) as hdf:
#        if run_string == 'onbeamR1':
#            evtdf = onbeamR1_evt_df_sb
#            hdr = pd.read_hdf(offbeamfile_Run1, "hdr")
#        if run_string == 'onbeamR2': 
#            evtdf = onbeamR2_evt_df_sb
#            hdr = pd.read_hdf(offbeamfile_Run2, "hdr")
#        try:
#            hdf.put(key='evt', value=evtdf, format="fixed")
#            hdf.put(key='hdr', value=hdr, format="fixed")
#        except Exception as e:
#            print("Table %s failed to save, skipping. Exception: %s" % ('hdr', str(e)))

#del offbeamR1_evt
#del offeamR1_evt_df
#del offbeamR1_evt_df_sb

#del offbeamR2_evt
#del offbeamR2_evt_df
#del offbeamR2_evt_df_sb

# ~~~~~~~ notes ~~~~~~~

#onbeamR1 livetime:  9712335.74
#onbeamR2 livetime:  40955958.70833331
#offbeamR1 livetime:  8327442.930000001 #this one went up! Hmm..
#offbeamR2 livetime:  39448777.45500001

#onbeamR1 POT:  4.651603646822754e+19
#onbeamR2 POT:  1.9447555727369336e+20

    
    # BEFORE CONCRETE FIX (GOOD RUNS LIST ALSO CHANGED THERE, CAUSING DIFFERENCE FOR DATA)
    
    # onbeamR1 livetime:   9882862.456666667
    # onbeamR2 livetime:  42174524.55000002
    # offbeamR1 livetime:  5028250.154999999
    # offbeamR2 livetime: 50284602.855000004

    # mccoh POT:  5.926336e+21
    #
    # onbeamR1 POT:  5.128076572325904e+19
    # onbeamR2 POT:  2.122704992632837e+20

#systFile:  /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_ind1bin0_evt.df
#mc_incoh livetime:  0.0
#mc_incoh POT:  4.821847e+21

#systFile:  /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_ind1bin14_evt.df
#mc_incoh livetime:  0.0
#mc_incoh POT:  5.5285587e+21

#systFile:  /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_ind0glo_evt.df
#mc_incoh livetime:  0.0
#mc_incoh POT:  4.1396564e+21

#systFile:  /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_ind0ghi_evt.df
#mc_incoh livetime:  0.0
#mc_incoh POT:  3.5396618e+21

#systFile:  /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_noiselhi_evt.df
#mc_incoh livetime:  0.0
#mc_incoh POT:  5.3094546e+21

#systFile:  /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_sce2x_evt.df
#mc_incoh livetime:  0.0
#mc_incoh POT:  5.1986385e+21

#systFile:  /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_ind0nom_evt.df
#mc_incoh livetime:  0.0
#mc_incoh POT:  5.3344124e+21

# ~~~~~~~ end of notes ~~~~~~~

        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# HERE IS THE STUFF TO RUN EVERY OTHER TIME:
    
    
    

for label in detVar_labels: print(label)

onbeamR1_pot = 4.651603646822754e+19
print('onbeamR1_pot: ', onbeamR1_pot)
onbeamR1_livetime = 9712335.74

onbeamR2_pot = 1.9447555727369336e+20
print('onbeamR2_pot: ', onbeamR2_pot)
onbeamR2_livetime = 40955958.70833331

onbeam_pot = onbeamR1_pot + onbeamR2_pot
print('Combined onbeam POT: ', onbeam_pot)
GOAL_POT = onbeam_pot

offbeamR1_livetime = 8327442.930000001
offbeamR2_livetime = 39448777.45500001

filedir_noConcrete = "/exp/icarus/data/users/jdyer/muon_datamc/samples_240909/lightened_before_ConcreteFix/" # this stuff restricted to sideband
filedir = "/exp/icarus/data/users/jdyer/muon_datamc/samples_240909/lightened_wConcreteFix_240909/" # this stuff restricted to sideband
all_mc_filedir = "/exp/icarus/data/users/jdyer/muon_datamc/samples_240909/allMC_ie_no_restriction_to_sideband_wConcreteFix_240909/" # this stuff NOT restricted to sideband.
# Note: The reason I am using these paths instead of savedir_sb and all_mc_filedir is becuase when reorganing this in Sept. 2024, I wanted to redo the above stuff for reproducability, and was planning to save it at the savedir_sb and all_mc_filedir. But, I couldn't complete that due to low disk space. The samples I had made previously were copied into filedir_noConcrete, filedir, and all_mc_filedir above though, so just keep using those. 


#~~~~~~~~~~~~~

# NOTE: need to run with simpler_add_hdr_info un-commented if I want to be able to look up relevant info needed to locate reco1 files for event displays.

print(h5py.File(filedir + "mc_incoh", 'r').keys())
#mc_incoh_evt = pd.read_hdf(mc_filedir + "mc_incoh", "evt")
mc_incoh_evt = pd.read_hdf(filedir + "mc_incoh", "evt")
#print(mc_inccoh_evt.shape)
sideband_mask = sb_numi_angle_mask(mc_incoh_evt)[0] & sb_open_angle_mask(mc_incoh_evt)[0]
#mc_incoh_evt = mc_incoh_evt[sideband_mask]
mc_incoh_evt["mc"] = True
mc_incoh_evt["mc_incoh"] = True
mc_incoh_evt["mccoh"] = False
mc_incoh_evt["onbeam"] = False
mc_incoh_evt["offbeam"] = False
mc_incoh_evt["Run1"] = False
mc_incoh_evt["Run2"] = False
#print(mc_inccoh_evt.shape)
mc_incoh_hdr = pd.read_hdf(filedir + "mc_incoh", "hdr")
mc_incoh_pot = np.sum(mc_incoh_hdr.pot * mc_incoh_hdr.first_in_subrun)
print('mc_incoh_pot: ', mc_incoh_pot)
mc_incoh_evt["scale"] = GOAL_POT/mc_incoh_pot
#mc_incoh_evt = simpler_add_hdr_info(mc_incoh_evt, mc_incoh_hdr)

print(h5py.File(filedir + "mccoh", 'r').keys())
#mccoh_evt = pd.read_hdf(mc_filedir + "mccoh", "evt")
mccoh_evt = pd.read_hdf(filedir + "mccoh", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_evt)[0] & sb_open_angle_mask(mccoh_evt)[0]
#mccoh_evt = mccoh_evt[sideband_mask]
mccoh_evt["mc"] = True
mccoh_evt["mc_incoh"] = False
mccoh_evt["mccoh"] = True
mccoh_evt["onbeam"] = False
mccoh_evt["offbeam"] = False
mccoh_evt["Run1"] = False
mccoh_evt["Run2"] = False
mccoh_hdr = pd.read_hdf(filedir + "mccoh", "hdr")
mccoh_hdr_pot = np.sum(mccoh_hdr.pot * mccoh_hdr.first_in_subrun)
print('mccoh_hdr_pot: ', mccoh_hdr_pot)
mccoh_evt["scale"] = GOAL_POT/mccoh_hdr_pot
#mccoh_evt = simpler_add_hdr_info(mccoh_evt, mccoh_hdr)

# ~~~~~~~~~~~~~~~~~~~
syst_evtdfs = []
syst_hdrs = []

#Middle Ind. Opaque
#Middle Ind. Transparent
#Front Ind. Gain Low
#Front Ind. Gain High
#Noise 1.2x
#Space Charge 2x
#Ind0 Nom

#print(h5py.File(filedir + "Middle_Ind._Opaque", 'r').keys())
mccoh_ind1bin0_evt = pd.read_hdf(filedir + "Middle_Ind._Opaque", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_ind1bin0_evt)[0] & sb_open_angle_mask(mccoh_ind1bin0_evt)[0]
mccoh_ind1bin0_evt = mccoh_ind1bin0_evt[sideband_mask]
mccoh_ind1bin0_evt["mc"] = True
mccoh_ind1bin0_evt["mc_incoh"] = False
mccoh_ind1bin0_evt["mccoh"] = True
mccoh_ind1bin0_evt["onbeam"] = False
mccoh_ind1bin0_evt["offbeam"] = False
mccoh_ind1bin0_evt["Run1"] = False
mccoh_ind1bin0_evt["Run2"] = False
mccoh_ind1bin0_evt["phi_NuMI_mcs"] = phi_NuMI(mccoh_ind1bin0_evt.trunk.trk, mccoh_ind1bin0_evt.branch.trk, BEAMDIR, 'mcs')[1]
syst_evtdfs.append(mccoh_ind1bin0_evt)
mccoh_ind1bin0_hdr = pd.read_hdf(filedir + "Middle_Ind._Opaque", "hdr")
mccoh_ind1bin0_hdr_pot = np.sum(mccoh_ind1bin0_hdr.pot * mccoh_ind1bin0_hdr.first_in_subrun)
syst_hdrs.append(mccoh_ind1bin0_hdr)
mccoh_ind1bin0_evt["scale"] = GOAL_POT/mccoh_ind1bin0_hdr_pot
#mccoh_ind1bin0_evt = simpler_add_hdr_info(mccoh_ind1bin0_evt, mccoh_ind1bin0_hdr)

#print(h5py.File(filedir + "Middle_Ind._Transparent", 'r').keys())
mccoh_ind1bin14_evt = pd.read_hdf(filedir + "Middle_Ind._Transparent", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_ind1bin14_evt)[0] & sb_open_angle_mask(mccoh_ind1bin14_evt)[0]
mccoh_ind1bin14_evt = mccoh_ind1bin14_evt[sideband_mask]
mccoh_ind1bin14_evt["mc"] = True
mccoh_ind1bin14_evt["mc_incoh"] = False
mccoh_ind1bin14_evt["mccoh"] = True
mccoh_ind1bin14_evt["onbeam"] = False
mccoh_ind1bin14_evt["offbeam"] = False
mccoh_ind1bin14_evt["Run1"] = False
mccoh_ind1bin14_evt["Run2"] = False
mccoh_ind1bin14_evt["phi_NuMI_mcs"] = phi_NuMI(mccoh_ind1bin14_evt.trunk.trk, mccoh_ind1bin14_evt.branch.trk, BEAMDIR, 'mcs')[1]
syst_evtdfs.append(mccoh_ind1bin14_evt)
mccoh_ind1bin14_hdr = pd.read_hdf(filedir + "Middle_Ind._Transparent", "hdr")
mccoh_ind1bin14_hdr_pot = np.sum(mccoh_ind1bin14_hdr.pot * mccoh_ind1bin14_hdr.first_in_subrun)
syst_hdrs.append(mccoh_ind1bin14_hdr)
mccoh_ind1bin14_evt["scale"] = GOAL_POT/mccoh_ind1bin14_hdr_pot
#mccoh_ind1bin14_evt = simpler_add_hdr_info(mccoh_ind1bin14_evt, mccoh_ind1bin14_hdr)

#print(h5py.File(filedir + "Front_Ind._Gain_Low", 'r').keys())
mccoh_ind0gainlo_evt = pd.read_hdf(filedir + "Front_Ind._Gain_Low", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_ind0gainlo_evt)[0] & sb_open_angle_mask(mccoh_ind0gainlo_evt)[0]
mccoh_ind0gainlo_evt = mccoh_ind0gainlo_evt[sideband_mask]
mccoh_ind0gainlo_evt["mc"] = True
mccoh_ind0gainlo_evt["mc_incoh"] = False
mccoh_ind0gainlo_evt["mccoh"] = True
mccoh_ind0gainlo_evt["onbeam"] = False
mccoh_ind0gainlo_evt["offbeam"] = False
mccoh_ind0gainlo_evt["Run1"] = False
mccoh_ind0gainlo_evt["Run2"] = False
mccoh_ind0gainlo_evt["phi_NuMI_mcs"] = phi_NuMI(mccoh_ind0gainlo_evt.trunk.trk, mccoh_ind0gainlo_evt.branch.trk, BEAMDIR, 'mcs')[1]
syst_evtdfs.append(mccoh_ind0gainlo_evt)
mccoh_ind0gainlo_hdr = pd.read_hdf(filedir + "Front_Ind._Gain_Low", "hdr")
mccoh_ind0gainlo_hdr_pot = np.sum(mccoh_ind0gainlo_hdr.pot * mccoh_ind0gainlo_hdr.first_in_subrun)
syst_hdrs.append(mccoh_ind0gainlo_hdr)
mccoh_ind0gainlo_evt["scale"] = GOAL_POT/mccoh_ind0gainlo_hdr_pot
#mccoh_ind0gainlo_evt = simpler_add_hdr_info(mccoh_ind0gainlo_evt, mccoh_ind0gainlo_hdr)

#print(h5py.File(filedir + "Front_Ind._Gain_High", 'r').keys())
mccoh_ind0gainhi_evt = pd.read_hdf(filedir + "Front_Ind._Gain_High", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_ind0gainhi_evt)[0] & sb_open_angle_mask(mccoh_ind0gainhi_evt)[0]
mccoh_ind0gainhi_evt = mccoh_ind0gainhi_evt[sideband_mask]
mccoh_ind0gainhi_evt["mc"] = True
mccoh_ind0gainhi_evt["mc_incoh"] = False
mccoh_ind0gainhi_evt["mccoh"] = True
mccoh_ind0gainhi_evt["onbeam"] = False
mccoh_ind0gainhi_evt["offbeam"] = False
mccoh_ind0gainhi_evt["Run1"] = False
mccoh_ind0gainhi_evt["Run2"] = False
mccoh_ind0gainhi_evt["phi_NuMI_mcs"] = phi_NuMI(mccoh_ind0gainhi_evt.trunk.trk, mccoh_ind0gainhi_evt.branch.trk, BEAMDIR, 'mcs')[1]
syst_evtdfs.append(mccoh_ind0gainhi_evt)
mccoh_ind0gainhi_hdr = pd.read_hdf(filedir + "Front_Ind._Gain_High", "hdr")
mccoh_ind0gainhi_hdr_pot = np.sum(mccoh_ind0gainhi_hdr.pot * mccoh_ind0gainhi_hdr.first_in_subrun)
syst_hdrs.append(mccoh_ind0gainhi_hdr)
mccoh_ind0gainhi_evt["scale"] = GOAL_POT/mccoh_ind0gainhi_hdr_pot
#mccoh_ind0gainhi_evt = simpler_add_hdr_info(mccoh_ind0gainhi_evt, mccoh_ind0gainhi_hdr)

#print(h5py.File(filedir + "Noise_1.2x", 'r').keys())
mccoh_noisehi_evt = pd.read_hdf(filedir + "Noise_1.2x", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_noisehi_evt)[0] & sb_open_angle_mask(mccoh_noisehi_evt)[0]
mccoh_noisehi_evt = mccoh_noisehi_evt[sideband_mask]
mccoh_noisehi_evt["mc"] = True
mccoh_noisehi_evt["mc_incoh"] = False
mccoh_noisehi_evt["mccoh"] = True
mccoh_noisehi_evt["onbeam"] = False
mccoh_noisehi_evt["offbeam"] = False
mccoh_noisehi_evt["Run1"] = False
mccoh_noisehi_evt["Run2"] = False
mccoh_noisehi_evt["phi_NuMI_mcs"] = phi_NuMI(mccoh_noisehi_evt.trunk.trk, mccoh_noisehi_evt.branch.trk, BEAMDIR, 'mcs')[1]
syst_evtdfs.append(mccoh_noisehi_evt)
mccoh_noisehi_hdr = pd.read_hdf(filedir + "Noise_1.2x", "hdr")
mccoh_noisehi_hdr_pot = np.sum(mccoh_noisehi_hdr.pot * mccoh_noisehi_hdr.first_in_subrun)
syst_hdrs.append(mccoh_noisehi_hdr)
mccoh_noisehi_evt["scale"] = GOAL_POT/mccoh_noisehi_hdr_pot
#mccoh_noisehi_evt = simpler_add_hdr_info(mccoh_noisehi_evt, mccoh_noisehi_hdr)

#print(h5py.File(filedir + "Space_Charge_2x", 'r').keys())
mccoh_sc2_evt = pd.read_hdf(filedir + "Space_Charge_2x", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_sc2_evt)[0] & sb_open_angle_mask(mccoh_sc2_evt)[0]
mccoh_sc2_evt = mccoh_sc2_evt[sideband_mask]
mccoh_sc2_evt["mc"] = True
mccoh_sc2_evt["mc_incoh"] = False
mccoh_sc2_evt["mccoh"] = True
mccoh_sc2_evt["onbeam"] = False
mccoh_sc2_evt["offbeam"] = False
mccoh_sc2_evt["Run1"] = False
mccoh_sc2_evt["Run2"] = False
mccoh_sc2_evt["phi_NuMI_mcs"] = phi_NuMI(mccoh_sc2_evt.trunk.trk, mccoh_sc2_evt.branch.trk, BEAMDIR, 'mcs')[1]
syst_evtdfs.append(mccoh_sc2_evt)
mccoh_sc2_hdr = pd.read_hdf(filedir + "Space_Charge_2x", "hdr")
mccoh_sc2_hdr_pot = np.sum(mccoh_sc2_hdr.pot * mccoh_sc2_hdr.first_in_subrun)
syst_hdrs.append(mccoh_sc2_hdr)
mccoh_sc2_evt["scale"] = GOAL_POT/mccoh_sc2_hdr_pot
#mccoh_sc2_evt = simpler_add_hdr_info(mccoh_sc2_evt, mccoh_sc2_hdr)

#print(h5py.File(filedir + "Ind0_Nom", 'r').keys())
mccoh_ind0nom_evt = pd.read_hdf(filedir + "Ind0_Nom", "evt")
sideband_mask = sb_numi_angle_mask(mccoh_ind0nom_evt)[0] & sb_open_angle_mask(mccoh_ind0nom_evt)[0]
mccoh_ind0nom_evt = mccoh_ind0nom_evt[sideband_mask]
mccoh_ind0nom_evt["mc"] = True
mccoh_ind0nom_evt["mc_incoh"] = False
mccoh_ind0nom_evt["mccoh"] = True
mccoh_ind0nom_evt["onbeam"] = False
mccoh_ind0nom_evt["offbeam"] = False
mccoh_ind0nom_evt["Run1"] = False
mccoh_ind0nom_evt["Run2"] = False
mccoh_ind0nom_evt["phi_NuMI_mcs"] = phi_NuMI(mccoh_ind0nom_evt.trunk.trk, mccoh_ind0nom_evt.branch.trk, BEAMDIR, 'mcs')[1]
syst_evtdfs.append(mccoh_ind0nom_evt)
mccoh_ind0nom_hdr = pd.read_hdf(filedir + "Ind0_Nom", "hdr")
mccoh_ind0nom_hdr_pot = np.sum(mccoh_ind0nom_hdr.pot * mccoh_ind0nom_hdr.first_in_subrun)
syst_hdrs.append(mccoh_ind0nom_hdr)
mccoh_ind0nom_evt["scale"] = GOAL_POT/mccoh_ind0nom_hdr_pot
#mccoh_ind0nom_evt = simpler_add_hdr_info(mccoh_ind0nom_evt, mccoh_ind0nom_hdr)

# Combine detector variation samples into a list:

new_syst_evtdfs = []
for df in syst_evtdfs:
    new_df = df[satisfies_new_FV(df)]
    new_syst_evtdfs.append(add_calculated_evtdf_cols(new_df))
syst_evtdfs = new_syst_evtdfs

# ~~~~~~~~~~~~~~~~~~~
# DATA

#print(h5py.File(filedir + "onbeamR1", 'r').keys())
onbeamR1_evt = pd.read_hdf(filedir + "onbeamR1", "evt")
sideband_mask = sb_numi_angle_mask(onbeamR1_evt)[0] & sb_open_angle_mask(onbeamR1_evt)[0]
onbeamR1_evt = onbeamR1_evt[sideband_mask]
onbeamR1_evt["mc"] = False
onbeamR1_evt["mc_incoh"] = False
onbeamR1_evt["mccoh"] = False
onbeamR1_evt["onbeam"] = True
onbeamR1_evt["offbeam"] = False
onbeamR1_evt["Run1"] = True
onbeamR1_evt["Run2"] = False
onbeamR1_hdr = pd.read_hdf(filedir + "onbeamR1", "hdr")
onbeamR1_evt["scale"] = 1
#onbeamR1_evt = simpler_add_hdr_info(onbeamR1_evt, onbeamR1_hdr)


#print(h5py.File(filedir + "onbeamR2", 'r').keys())
onbeamR2_evt = pd.read_hdf(filedir + "onbeamR2", "evt")
sideband_mask = sb_numi_angle_mask(onbeamR2_evt)[0] & sb_open_angle_mask(onbeamR2_evt)[0]
onbeamR2_evt = onbeamR2_evt[sideband_mask]
onbeamR2_evt["mc"] = False
onbeamR2_evt["mc_incoh"] = False
onbeamR2_evt["mccoh"] = False
onbeamR2_evt["onbeam"] = True
onbeamR2_evt["offbeam"] = False
onbeamR2_evt["Run1"] = False
onbeamR2_evt["Run2"] = True
onbeamR2_hdr = pd.read_hdf(filedir + "onbeamR2", "hdr")
onbeamR2_evt["scale"] = 1
#onbeamR2_evt = simpler_add_hdr_info(onbeamR2_evt, onbeamR2_hdr)

#print(h5py.File(filedir + "offbeamR1", 'r').keys())
offbeamR1_evt = pd.read_hdf(filedir + "offbeamR1", "evt")
sideband_mask = sb_numi_angle_mask(offbeamR1_evt)[0] & sb_open_angle_mask(offbeamR1_evt)[0]
offbeamR1_evt = offbeamR1_evt[sideband_mask]
offbeamR1_evt["mc"] = False
offbeamR1_evt["mc_incoh"] = False
offbeamR1_evt["mccoh"] = False
offbeamR1_evt["onbeam"] = False
offbeamR1_evt["offbeam"] = True
offbeamR1_evt["Run1"] = True
offbeamR1_evt["Run2"] = False
offbeamR1_hdr = pd.read_hdf(filedir + "offbeamR1", "hdr")
offbeamR1_evt["scale"] = onbeamR1_livetime/offbeamR1_livetime
#offbeamR1_evt = simpler_add_hdr_info(offbeamR1_evt, offbeamR1_hdr)

#print(h5py.File(filedir + "offbeamR2", 'r').keys())
offbeamR2_evt = pd.read_hdf(filedir + "offbeamR2", "evt")
sideband_mask = sb_numi_angle_mask(offbeamR2_evt)[0] & sb_open_angle_mask(offbeamR2_evt)[0]
offbeamR2_evt = offbeamR2_evt[sideband_mask]
offbeamR2_evt["mc"] = False
offbeamR2_evt["mc_incoh"] = False
offbeamR2_evt["mccoh"] = False
offbeamR2_evt["onbeam"] = False
offbeamR2_evt["offbeam"] = True
offbeamR2_evt["Run1"] = False
offbeamR2_evt["Run2"] = True
offbeamR2_hdr = pd.read_hdf(filedir + "offbeamR2", "hdr")
offbeamR2_evt["scale"] = onbeamR2_livetime/offbeamR2_livetime
#offbeamR2_evt = simpler_add_hdr_info(offbeamR2_evt, offbeamR2_hdr)































