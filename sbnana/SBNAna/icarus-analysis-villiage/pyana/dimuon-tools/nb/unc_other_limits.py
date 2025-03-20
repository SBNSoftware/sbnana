# Load in the python data objects for the other limits, so they're ready to plot in my sensitivity curves.

import numpy as np

### EXISTING HPS LIMITS:

other_limits_dir = '/exp/icarus/data/users/jdyer/muon_selection/other_limits/hps/'
# ICARUS Contained dimuon search
hps_Gray_x = [] # MeV
hps_Gray_y = [] # theta^2
with open(other_limits_dir+'data_from_plotDigitizer/Grays_contained_HPS_90CL_data', 'r') as file:
    for line in file:
        x, y = line.strip().split()
        hps_Gray_x.append(float(x))
        hps_Gray_y.append(float(y))    
# NA62 search
hps_na62_x_lower_limits = [] # GeV
hps_na62_y_lower_limits = [] # sin^2(theta)
with open(other_limits_dir+'NA62_90CL_data_fromSource_lower_limits', 'r') as file:
    for line in file:
        x, y = line.strip().split()
        hps_na62_x_lower_limits.append(float(x)*1000.) # MeV
        hps_na62_y_lower_limits.append(float(y))
hps_na62_y_lower_limits = list(np.arcsin(np.sqrt(hps_na62_y_lower_limits))**2) # theta^2

hps_na62_x_upper_limits = [] # GeV
hps_na62_y_upper_limits = [] # sin^2(theta)
with open(other_limits_dir+'NA62_90CL_data_fromSource_upper_limits', 'r') as file:
    for line in file:
        x, y = line.strip().split()
        hps_na62_x_upper_limits.append(float(x)*1000.) # MeV
        hps_na62_y_upper_limits.append(float(y))
hps_na62_x_upper_limits.reverse()
hps_na62_y_upper_limits.reverse()
hps_na62_y_upper_limits = list(np.arcsin(np.sqrt(hps_na62_y_upper_limits))**2) # theta^2

hps_na62_x = hps_na62_x_lower_limits + hps_na62_x_upper_limits
hps_na62_y = hps_na62_y_lower_limits + hps_na62_y_upper_limits

# MicroBooNE search
hps_uB_x = []
hps_uB_y = []
with open(other_limits_dir+'uBSensitivity.txt', 'r') as file:
    for line in file:
        try:
            x, y = line.strip().split(', ')
            hps_uB_x.append(float(x))
            hps_uB_y.append(float(y))
        except:
            null_var = 0
            #print('this line doesnt look like data: \n', line)
# E949 Search
hps_E949_x = []
hps_E949_y = []
with open(other_limits_dir+'E949Sensitivity.txt', 'r') as file:
    for line in file:
        try:
            x, y = line.strip().split(', ')
            hps_E949_x.append(float(x))
            hps_E949_y.append(float(y))
        except:
            null_var = 0
            #print('This line doesnt look like data: \n', line)
# LHCb Search
hps_LHCb_x = []
hps_LHCb_y = []
with open(other_limits_dir+'LHCbSensitivity.txt', 'r') as file:
    for line in file:
        try:
            x, y = line.strip().split(', ')
            hps_LHCb_x.append(float(x))
            hps_LHCb_y.append(float(y))
        except:
            null_var = 0
            #print('this line doesnt look like data: \n', line)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
### EXISTING HPS LIMITS:

# convention: y = 1/fa
# save masses in GeV

# -----

# CODOMINANCE, RUNNING COUPLING for c_mu

otherLims_runningCmu_dir = '/exp/icarus/data/users/jdyer/muon_selection/other_limits/alps/'

# ICARUS Contained search
masses_lowlim = [] # save it as MeV
inv_fa_lowlim = []
with open(otherLims_runningCmu_dir+'Grays_obsLimit_CD_uurun_lowfa', 'r') as file:
    for line in file:
        a, b = line.strip().split(', ')
        masses_lowlim.append(float(a))
        inv_fa_lowlim.append(float(b))
masses_hilim = [] # save it as MeV
inv_fa_hilim = []
with open(otherLims_runningCmu_dir+'Grays_obsLimit_CD_uurun_hifa', 'r') as file:
    for line in file:
        a, b = line.strip().split(', ')
        masses_hilim.append(float(a))
        inv_fa_hilim.append(float(b))
alp_Gray_x = np.concatenate(( np.array(masses_lowlim), np.flip(np.array(masses_hilim)) ))
alp_Gray_y = np.concatenate(( np.array(inv_fa_lowlim), np.flip(np.array(inv_fa_hilim)) ))

# CHARM gg - from Kaons
masses = [] # save it as MeV
fa_max = []
fa_min = []
with open(otherLims_runningCmu_dir+'limits_CHARM_CHARM_gg_CD_N4.txt', 'r') as file:
    for line in file:
        if '#' in line: 
            continue
        else:
            a, b, c = line.strip().split('\t')
            masses.append(float(a))
            fa_max.append(float(b))
            fa_min.append(float(c))
alp_charm_gg_x = np.concatenate(( np.array(masses), np.flip(np.array(masses)) ))
alp_charm_gg_y = 1./np.concatenate(( np.array(fa_max), np.flip(np.array(fa_min)) )) # 1/fa
alp_charm_gg_x = alp_charm_gg_x[np.isfinite(alp_charm_gg_y)]
alp_charm_gg_y = alp_charm_gg_y[np.isfinite(alp_charm_gg_y)]

# CHARM uu - from Kaons
masses = [] # save it as MeV
fa_max = []
fa_min = []
with open(otherLims_runningCmu_dir+'limits_CHARM_CHARM_uurun_CD_N2.txt', 'r') as file:
    for line in file:
        if '#' in line: 
            continue
        else:
            a, b, c = line.strip().split('\t')
            masses.append(float(a))
            fa_max.append(float(c))
            fa_min.append(float(b))
alp_charm_uu_x = np.concatenate(( np.array(masses), np.flip(np.array(masses)) ))
alp_charm_uu_y = 1./np.concatenate(( np.array(fa_max), np.flip(np.array(fa_min)) )) # 1/fa
alp_charm_uu_x = alp_charm_uu_x[np.isfinite(alp_charm_uu_y)]
alp_charm_uu_y = alp_charm_uu_y[np.isfinite(alp_charm_uu_y)]

# CHARM gg - from psuedoscalars
masses = [] # save it as MeV
fa_max = []
fa_min = []
with open(otherLims_runningCmu_dir+'limits_CHARM_CHARMP0_gg_CD_N4.txt', 'r') as file:
    for line in file:
        if '#' in line: 
            continue
        else:
            a, b, c = line.strip().split('\t')
            masses.append(float(a))
            fa_max.append(float(b))
            fa_min.append(float(c))
alp_charmP0_gg_x = np.concatenate(( np.array(masses), np.flip(np.array(masses)) ))
alp_charmP0_gg_y = 1./np.concatenate(( np.array(fa_max), np.flip(np.array(fa_min)) )) # 1/fa
alp_charmP0_gg_x = alp_charmP0_gg_x[np.isfinite(alp_charmP0_gg_y)]
alp_charmP0_gg_y = alp_charmP0_gg_y[np.isfinite(alp_charmP0_gg_y)]
# enable not connecting the upper and lower limits:
charmP0_gg_x_1 = np.array(masses)
charmP0_gg_y_1 = 1./np.array(fa_max)
charmP0_gg_x_1 = charmP0_gg_x_1[np.isfinite(charmP0_gg_y_1)]
charmP0_gg_y_1 = charmP0_gg_y_1[np.isfinite(charmP0_gg_y_1)]
charmP0_gg_x_2 = np.array(masses)
charmP0_gg_y_2 = 1./np.array(fa_min)
charmP0_gg_x_2 = charmP0_gg_x_2[np.isfinite(charmP0_gg_y_2)]
charmP0_gg_y_2 = charmP0_gg_y_2[np.isfinite(charmP0_gg_y_2)]


# CHARM uu - from psuedoscalars
masses = [] # save it as MeV
fa_max = []
fa_min = []
with open(otherLims_runningCmu_dir+'limits_CHARM_CHARMP0_uurun_CD_N2.txt', 'r') as file:
    for line in file:
        if '#' in line: 
            continue
        else:
            a, b, c = line.strip().split('\t')
            masses.append(float(a))
            fa_max.append(float(c))
            fa_min.append(float(b))
alp_charmP0_uu_x = np.concatenate(( np.array(masses), np.flip(np.array(masses)) ))
alp_charmP0_uu_y = 1./np.concatenate(( np.array(fa_max), np.flip(np.array(fa_min)) )) # 1/fa
alp_charmP0_uu_x = alp_charmP0_uu_x[np.isfinite(alp_charmP0_uu_y)]
alp_charmP0_uu_y = alp_charmP0_uu_y[np.isfinite(alp_charmP0_uu_y)]
# enable not connecting the upper and lower limits:
charmP0_uu_x_1 = np.array(masses)
charmP0_uu_y_1 = 1./np.array(fa_max)
charmP0_uu_x_1 = charmP0_uu_x_1[np.isfinite(charmP0_uu_y_1)]
charmP0_uu_y_1 = charmP0_uu_y_1[np.isfinite(charmP0_uu_y_1)]
charmP0_uu_x_2 = np.array(masses)
charmP0_uu_y_2 = 1./np.array(fa_min)
charmP0_uu_x_2 = charmP0_uu_x_2[np.isfinite(charmP0_uu_y_2)]
charmP0_uu_y_2 = charmP0_uu_y_2[np.isfinite(charmP0_uu_y_2)]

# NA62
# Red in Figure 9 right panel here: https://arxiv.org/pdf/2103.15389
# This is an alp-interpretation of same data for HPS.
masses = [] # save it as MeV
inv_fa_max = []
inv_fa_min = []
with open(otherLims_runningCmu_dir+'NA62_bound_CD_cuu0_mhi.txt', 'r') as file:
    for line in file:
        a, b, c = line.strip().split('\t')
        masses.append(float(a))
        inv_fa_max.append(float(c))
        inv_fa_min.append(float(b))
alp_NA62_x = np.concatenate(( np.array(masses), np.flip(np.array(masses)) ))
alp_NA62_y = np.concatenate(( np.array(inv_fa_min), np.flip(np.array(inv_fa_max)) )) # 1/fa

# uB
masses = [] # save it as MeV
fa_max = []
fa_min = []
with open(otherLims_runningCmu_dir+'uBsensitivity_ALP_CD.txt', 'r') as file:
    for line in file:
        if '#' in line: 
            continue
        else:
            a, b, c = line.strip().split('\t')
            masses.append(float(a))
            fa_max.append(float(b))
            fa_min.append(float(c))
alp_uB_x = np.concatenate(( np.array(masses), np.flip(np.array(masses)) ))
alp_uB_y = 1./np.concatenate(( np.array(fa_max), np.flip(np.array(fa_min)) )) # 1/fa














