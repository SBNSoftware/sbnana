import uproot
import pandas as pd 
import numpy as np 
import h5py
import math
from unc_funcs import *

## Note July 28, 2025: I've just reprocessed dataframes using the CAFs I have available on scratch to include more stuff about reconstructed daughters and the MCS scattering angles, as well as updated flux CV and systematic weights. Updated dataframes are missing for the/a nominal neutrino sample and all of the detector variation samples so far.

updated_df_dir = "/exp/icarus/data/users/jdyer/dimuon-data/2507_remake/"

# DATAFRAME FILES:

gray_df_dir = "/exp/icarus/data/users/gputnam/thesis-work/DMCP2023G/mc-F/"
detVar_labels = [
    ##"Middle Ind. Opaque",
    ##"Middle Ind. Transparent",
    "Middle Ind. Opacity",
    ##"Front Ind. Gain Low",
    ##"Front Ind. Gain High",
    #"Front Ind. Gain", # Uncomment later once I have it.
    "Noise 1.2x", 
    "Space Charge 2x",
    "Ind0 Nom" # ??
]

# Higgs Files
# Note: Gray made these. They exclude final state pions.
higgs_files = [
    #updated_df_dir + "higgs_250915/ms300.df",
    #updated_df_dir + "higgs_250915/ms320.df",
    #updated_df_dir + "higgs_250915/ms350.df"
    updated_df_dir + "ms220.df",
    updated_df_dir + "ms240.df",
    updated_df_dir + "ms260.df",
    updated_df_dir + "ms280.df",
    updated_df_dir + "ms300.df",
    updated_df_dir + "ms320.df",
    updated_df_dir + "ms340.df"
]

higgs_220_detVar_files = [
    gray_df_dir + "F2-Higgs_M220_ind1bin0_evt.df", # "Middle Ind. Opaque"
    gray_df_dir + "F2-Higgs_M220_ind1bin14_evt.df", # "Middle Ind. Transparent"
    gray_df_dir + "F2-Higgs_M220_ind0glo_evt.df", # "Front Ind. Gain Low"
    gray_df_dir + "F2-Higgs_M220_ind0ghi_evt.df", # "Front Ind. Gain High"
    gray_df_dir + "F2-Higgs_M220_noiselhi_evt.df", # "Noise 1.2x"
    gray_df_dir + "F2-Higgs_M220_sce2x_evt.df", # "Space Charge 2x"
    gray_df_dir + "F2-Higgs_M220_ind0nom_evt.df" # "Ind 0 Nom"
]
higgs_240_detVar_files = [
    gray_df_dir + "F2-Higgs_M240_ind1bin0_evt.df", # "Middle Ind. Opaque"
    gray_df_dir + "F2-Higgs_M240_ind1bin14_evt.df", # "Middle Ind. Transparent"
    gray_df_dir + "F2-Higgs_M240_ind0glo_evt.df", # "Front Ind. Gain Low"
    gray_df_dir + "F2-Higgs_M240_ind0ghi_evt.df", # "Front Ind. Gain High"
    gray_df_dir + "F2-Higgs_M240_noiselhi_evt.df", # "Noise 1.2x"
    gray_df_dir + "F2-Higgs_M240_sce2x_evt.df", # "Space Charge 2x"
    gray_df_dir + "F2-Higgs_M240_ind0nom_evt.df" # "Ind 0 Nom"
]
higgs_260_detVar_files = [
    gray_df_dir + "F2-Higgs_M260_ind1bin0_evt.df", # "Middle Ind. Opaque"
    gray_df_dir + "F2-Higgs_M260_ind1bin14_evt.df", # "Middle Ind. Transparent"
    gray_df_dir + "F2-Higgs_M260_ind0glo_evt.df", # "Front Ind. Gain Low"
    gray_df_dir + "F2-Higgs_M260_ind0ghi_evt.df", # "Front Ind. Gain High"
    gray_df_dir + "F2-Higgs_M260_noiselhi_evt.df", # "Noise 1.2x"
    gray_df_dir + "F2-Higgs_M260_sce2x_evt.df", # "Space Charge 2x"
    gray_df_dir + "F2-Higgs_M260_ind0nom_evt.df" # "Ind 0 Nom"
]
higgs_280_detVar_files = [
    gray_df_dir + "F2-Higgs_M280_ind1bin0_evt.df", # "Middle Ind. Opaque"
    gray_df_dir + "F2-Higgs_M280_ind1bin14_evt.df", # "Middle Ind. Transparent"
    gray_df_dir + "F2-Higgs_M280_ind0glo_evt.df", # "Front Ind. Gain Low"
    gray_df_dir + "F2-Higgs_M280_ind0ghi_evt.df", # "Front Ind. Gain High"
    gray_df_dir + "F2-Higgs_M280_noiselhi_evt.df", # "Noise 1.2x"
    gray_df_dir + "F2-Higgs_M280_sce2x_evt.df", # "Space Charge 2x"
    gray_df_dir + "F2-Higgs_M280_ind0nom_evt.df" # "Ind 0 Nom"
]
higgs_300_detVar_files = [
    gray_df_dir + "F2-Higgs_M300_ind1bin0_evt.df", # "Middle Ind. Opaque"
    gray_df_dir + "F2-Higgs_M300_ind1bin14_evt.df", # "Middle Ind. Transparent"
    gray_df_dir + "F2-Higgs_M300_ind0glo_evt.df", # "Front Ind. Gain Low"
    gray_df_dir + "F2-Higgs_M300_ind0ghi_evt.df", # "Front Ind. Gain High"
    gray_df_dir + "F2-Higgs_M300_noiselhi_evt.df", # "Noise 1.2x"
    gray_df_dir + "F2-Higgs_M300_sce2x_evt.df", # "Space Charge 2x"
    gray_df_dir + "F2-Higgs_M300_ind0nom_evt.df" # "Ind 0 Nom"
]
higgs_340_detVar_files = [
    gray_df_dir + "F2-Higgs_M340_ind1bin0_evt.df", # "Middle Ind. Opaque"
    gray_df_dir + "F2-Higgs_M340_ind1bin14_evt.df", # "Middle Ind. Transparent"
    gray_df_dir + "F2-Higgs_M340_ind0glo_evt.df", # "Front Ind. Gain Low"
    gray_df_dir + "F2-Higgs_M340_ind0ghi_evt.df", # "Front Ind. Gain High"
    gray_df_dir + "F2-Higgs_M340_noiselhi_evt.df", # "Noise 1.2x"
    gray_df_dir + "F2-Higgs_M340_sce2x_evt.df", # "Space Charge 2x"
    gray_df_dir + "F2-Higgs_M340_ind0nom_evt.df" # "Ind 0 Nom"
]

# ALP Files (ALPs WITHOUT MASS SUPPRESSION)

alp_nosup_files = [
    updated_df_dir + "ma300_faE6.df",
    #updated_df_dir + "ma340_faE6.df",
    updated_df_dir + "ma350_faE6.df",
    updated_df_dir + "ma400_faE6.df",
    updated_df_dir + "ma450_faE6.df",
    updated_df_dir + "ma500_faE6.df",
    updated_df_dir + "ma540_faE6.df",
    updated_df_dir + "ma600_faE6.df",
    updated_df_dir + "ma650_faE6.df",
    updated_df_dir + "ma700_faE6.df"
]

alp_300_detVar_files = [
    "/exp/icarus/data/users/jdyer/dimuon-data/2504_AlpDetectorVariations/Alp_M300_ind1bin0.df", # "Middle Ind. Opaque"
    "/exp/icarus/data/users/jdyer/dimuon-data/2504_AlpDetectorVariations/Alp_M300_ind1bin14.df", # "Middle Ind. Transparent"
    "/exp/icarus/data/users/jdyer/dimuon-data/2504_AlpDetectorVariations/Alp_M300_ind0glo.df", # "Front Ind. Gain Low"
    "/exp/icarus/data/users/jdyer/dimuon-data/2504_AlpDetectorVariations/Alp_M300_ind0ghi.df", # "Front Ind. Gain High"
    "/exp/icarus/data/users/jdyer/dimuon-data/2504_AlpDetectorVariations/Alp_M300_noiselhi.df", # "Noise 1.2x"
    "/exp/icarus/data/users/jdyer/dimuon-data/2504_AlpDetectorVariations/Alp_M300_sce2x.df", # "Space Charge 2x"
    "/exp/icarus/data/users/jdyer/dimuon-data/2504_AlpDetectorVariations/Alp_M300_ind0nom.df" # "Ind 0 Nom"
]

# Nu Files

nu_file = "/exp/icarus/data/users/jdyer/dimuon-data/2507_remake/nomNus.df" # This dataframe was made August 18, 2025. To make it, I staged all the files under the samweb definition "gputnam_DMCP2023F_F-MCNuPhase2_flatcafs," and ran the dataframe-making tool over all of those to make the dataframe. Many "nan index" warnings were thrown in the making of the dataframe, and it got made quickly. The resulting dataframe is also smaller than expected (only 2379 events), when the samweb definition file has [~6k] files associated with it. So, I definitely think that the dataframe doesn't have all the events that SHOULD be in the samweb definition (maybe some files were deprecated?). Anyways, this shouldn't cause any issues in principle, so use for now. It should help with stats.
#nu_file = "/exp/icarus/data/users/jdyer/dimuon-data/2507_remake/nomNus_251007spare.df"

#nu_file = gray_df_dir + "F-MCNuPhase2_evt.df"
#cohlike_nu_file = gray_df_dir + "F-CohLike_nom_evt.df"
#cohlike_nu_file_extra_2501 = "/exp/icarus/data/users/jdyer/dimuon-data/2501_moreCohlike/more_Cohlike.df" # reco1 files from these have been deleted from scratch.
#cohlike_nu_file_extra_2502_12345 = "/exp/icarus/data/users/jdyer/dimuon-data/2502_moreCohlike/12345/even_more_Cohlike_12345.df"
#cohlike_nu_file_extra_2502_12 = "/exp/icarus/data/users/jdyer/dimuon-data/2502_moreCohlike/12/even_more_Cohlike_12.df"
#cohlike_nu_file_extra_2502_34 = "/exp/icarus/data/users/jdyer/dimuon-data/2502_moreCohlike/34/even_more_Cohlike_34.df"
#cohlike_nu_file_extra_2502_5 = "/exp/icarus/data/users/jdyer/dimuon-data/2502_moreCohlike/5/even_more_Cohlike_5.df"
#cohlike_nu_file_extra_2502_6789 = "/exp/icarus/data/users/jdyer/dimuon-data/2502_moreCohlike/6789/even_more_Cohlike_6789.df"

# note: This sample contains coh-LIKE events. Really, it is any event with a muon or a pion in the final state. This means that the interaction isn't necessarily a coherent interaction, but it just looks like one.

# note: It is possible that there is pileup in the neutrino interactions. If it happens that two neutrino interactions overlap in a single reconstructed slice, then even if the "bonus" one is not coh-like, it would be included in this sample as well. Having two neutrino interactions overlap is somewhat likely in NuMI.

#cohlike_nu_file_extra_2501 = updated_df_dir + "more_Cohlike.df"
#cohlike_nu_file_extra_2502_12 = updated_df_dir + "even_more_Cohlike_12.df"
#cohlike_nu_file_extra_2502_5 = updated_df_dir + "even_more_Cohlike_5.df"
#cohlike_nu_file_extra_2502_6789 = updated_df_dir + "even_more_Cohlike_6789.df"

cohlike_1 = updated_df_dir+"cohlike-F-1.df"
cohlike_2 = updated_df_dir+"cohlike-F-2.df"

## note: We only have detector variation samples for coherent-like events. 
#cohlike_nu_detVar_files = [
#    gray_df_dir + "F-CohLike_ind1bin0_evt.df", # "Middle Ind. Opaque"
#    gray_df_dir + "F-CohLike_ind1bin14_evt.df", # "Middle Ind. Transparent"
#    gray_df_dir + "F-CohLike_ind0glo_evt.df", # "Front Ind. Gain Low"
#    gray_df_dir + "F-CohLike_ind0ghi_evt.df", # "Front Ind. Gain High"
#    gray_df_dir + "F-CohLike_noiselhi_evt.df", # "Noise 1.2x"
#    gray_df_dir + "F-CohLike_sce2x_evt.df", # "Space Charge 2x"
#    gray_df_dir + "F-CohLike_ind0nom_evt.df" # "Ind 0 Nom"
#    # MISSING: CALHI, CALLO # MAY 3, 2024 # May 5: Nevermind! Gray made that spectrum for events in their signal box.
#    ]

cohlike_nu_detVar_files = [
    updated_df_dir + "cohlike_ind1bin0-1.df", updated_df_dir + "cohlike_ind1bin0-2.df", # "Middle Ind. Opaque"
    #updated_df_dir + "cohlike_ind1bin14.df", # "Middle Ind. Transparent"
    #updated_df_dir + "cohlike_ind0glo.df", # "Front Ind. Gain Low"
    #updated_df_dir + "cohlike_ind0ghi.df", # "Front Ind. Gain Hi"
    updated_df_dir + "cohlike_noiselhi-1.df", updated_df_dir + "cohlike_noiselhi-2.df", # "Noise 1.2x"
    updated_df_dir + "cohlike_sce2x-1.df", updated_df_dir + "cohlike_sce2x-2.df", # "Space Charge 2x"
    updated_df_dir + "cohlike_ind0nom-1.df", updated_df_dir + "cohlike_ind0nom-2.df" # "Ind 0 Nom"

                          ]


# Data Files

onbeamfile_Run1 = gray_df_dir + "F-Run1-OnBeam_evt.df"
onbeamfile_Run2 = gray_df_dir + "F4-Run2-OnBeam_evt.df"

offbeamfile_Run1 = gray_df_dir + "F-Run1-OffBeam_evt.df"
offbeamfile_Run2 = gray_df_dir + "F2-Run2-OffBeam_evt.df"

# Notes:
# There is overlap between files that are used in different notebooks, but the data can be used differently. E.g.:
# For MC Event Selection and MC Studies, event data frames are extracted w/ Pandas function:
#    pd.read_hdf(filename, key='evt') 
# For Data/MC Comparisons nb, event data frames are extracted with custom function: 
#    data.mc_dataset(mccohfile, "evt", mcnukey="mcnuweight").df
# This function takes a while, becuase it uses information in [] to compute weights needed to run systematics, and then inserts that into the returned dataframe. 

