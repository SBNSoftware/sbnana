import numpy as np

# Load in the relative ALP fluxes Josh made w/ pythia, g4 primaries only, and g4 w secondaries:

rate_file = "/exp/icarus/data/users/jdyer/muon_selection/alp_flux_uncertainty/rate_table_logbin.txt"
#rate_file = "/exp/icarus/data/users/jdyer/muon_selection/alp_flux_uncertainty/rate_table_pctbin.txt"
#rate_file = "/exp/icarus/data/users/jdyer/muon_selection/alp_flux_uncertainty/rate_table.txt"
file_mass = []
E_low_bin_bound = [] 
E_hi_bin_bound = [] 
R_pythia = [] # only considers primaries
dR_pythia = []
R_g4_primOnly = []
dR_g4_primOnly = []
R_g4_wSec = []
dR_g4_wSec = []
not_first_line = False
with open(rate_file, 'r') as file:
    for line in file:
        if not_first_line:
            this_line = line.strip().split('\t')
            file_mass.append(int(this_line[0]))
            E_low_bin_bound.append(float(this_line[1])) 
            E_hi_bin_bound.append(float(this_line[2])) 
            R_pythia.append(float(this_line[3])) 
            dR_pythia.append(float(this_line[4])) 
            R_g4_primOnly.append(float(this_line[5])) 
            dR_g4_primOnly.append(float(this_line[6])) 
            R_g4_wSec.append(float(this_line[7])) 
            dR_g4_wSec.append(float(this_line[8]))
        else:
            not_first_line = True
E_low_bin_bound = np.array(E_low_bin_bound)
E_hi_bin_bound = np.array(E_hi_bin_bound)
R_pythia = np.array(R_pythia) 
R_g4_primOnly = np.array(R_g4_primOnly) 
R_g4_wSec = np.array(R_g4_wSec) 

# separate out by ALP mass:

E_low_bin_bounds = []
E_hi_bin_bounds = []
R_pythias = []
R_g4_primOnlys = []
R_g4_wSecs = []
for um in np.unique(np.array(file_mass)):
    m_mask = (np.array(file_mass)==um)
    E_low_bin_bounds.append(E_low_bin_bound[m_mask])
    E_hi_bin_bounds.append(E_hi_bin_bound[m_mask])
    R_pythias.append(R_pythia[m_mask])
    R_g4_primOnlys.append(R_g4_primOnly[m_mask])
    R_g4_wSecs.append(R_g4_wSec[m_mask])

test_dict = dict(zip(np.unique(file_mass)/1000., range(len(np.unique(file_mass)))))
    
def ALP_flux_rw_for_secondaries(alpE, alpM): # This is for a CV reweighting based on g4_primaries_only/g4_w_secondaries
    alp_mass_ind = test_dict[alpM]
    #print('alpE: ', alpE)
    E_bin_bounds = np.concatenate((E_low_bin_bounds[alp_mass_ind], np.array([E_hi_bin_bounds[alp_mass_ind][-1]])))
    try:
        lower_ind = np.where(E_bin_bounds>alpE)[0][0]-1
        upper_ind = np.where(E_bin_bounds>alpE)[0][0]
    except IndexError: # the ALP energy is higher than the max of the range we have rates for.
        # In this case, treat last energy bin as overflow bin
        lower_ind = len(E_bin_bounds)-2
        upper_ind = len(E_bin_bounds)-1
    Ng4_prim = R_g4_primOnlys[alp_mass_ind][lower_ind]
    Ng4_sec = R_g4_wSecs[alp_mass_ind][lower_ind]
    r = Ng4_sec/Ng4_prim
    #print(alpE, w, sep=', ')
    return(r)

def get_alp_flux_w(alpE, alp_mass_ind): 
# this is to get an uncertainty on the ALP flux based on g4_primaries_only/pythia_primaries_only
# these are reweights based on the ratio, so when applied to the CV-reweighted alps reweighted w/ above func, they will still produce the correct fractional uncertainty, which based on how the code works will also churn out the correct overall uncertainty.
    #print('alpE: ', alpE)
    E_bin_bounds = np.concatenate((E_low_bin_bounds[alp_mass_ind], np.array([E_hi_bin_bounds[alp_mass_ind][-1]])))
    try:
        lower_ind = np.where(E_bin_bounds>alpE)[0][0]-1
        upper_ind = np.where(E_bin_bounds>alpE)[0][0]
    except IndexError: # the ALP energy is higher than the max of the range we have rates for.
        # In this case, treat last energy bin as overflow bin
        lower_ind = len(E_bin_bounds)-2
        upper_ind = len(E_bin_bounds)-1
    Npythia = R_pythias[alp_mass_ind][lower_ind]
    Ng4 = R_g4_primOnlys[alp_mass_ind][lower_ind]
    w = Ng4/Npythia - 1
    #print(alpE, w, sep=', ')
    return(w)