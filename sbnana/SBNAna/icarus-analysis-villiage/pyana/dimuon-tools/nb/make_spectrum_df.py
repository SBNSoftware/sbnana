import sys
import os
import argparse
import uproot
import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np 
import matplotlib as mpl
import timewgt

from util import *
import var
import cut
import data
import hist

import importlib

from pyanalib import panda_helpers

parser = argparse.ArgumentParser(
    prog="make_spectrum_df",
    description="Convert an evt dataframe to a spectrum dataframe.",
)
parser.add_argument("filename")
parser.add_argument("-gl", "--gainlo", action="store_true", help="Gain Low Variation")
parser.add_argument("-gh", "--gainhi", action="store_true", help="Gain High Variation")
parser.add_argument("-cl", "--callo", action="store_true", help="dE/dx Low Variation")
parser.add_argument("-ch", "--calhi", action="store_true", help="dE/dx High Variation")
parser.add_argument("-ml", "--mcslo", action="store_true", help="MCS Low Variation")
parser.add_argument("-mh", "--mcshi", action="store_true", help="MCS High Variation")

parser.add_argument("-l", "--loose", action="store_true", help="Loose thumi Cut")
parser.add_argument("-rc", "--remove_cohlike", action="store_true", help="Remove coh-like events.")
parser.add_argument("-csb", "--cut_signal_box", action="store_true", help="Remove events inside signal box region.")
parser.add_argument("-cst", "--correct_split_tracks", action="store_true", help="Apply track splitting corrections w/ uncertainty.")
parser.add_argument("-d", "--data", action="store_true", help="File is data.")
parser.add_argument("-o", "--overwrite", action="store_true", help="Whether to overwrite existing spectra file.")
parser.add_argument("-od", "--output-dir", help="Output directory location")

NU_PER_POT_BASE = 1.0448039251186124e-16
parser.add_argument("-fp", "--fix-pot", help="Whether to fix POT to Coh-like amount (%f)." % NU_PER_POT_BASE, action="store_true")

args = parser.parse_args()

# get the type of any variation
if args.gainlo:
    vartype = "gainlo_"
elif args.gainhi:
    vartype = "gainhi_"
elif args.callo:
    vartype = "callo_"
elif args.calhi:
    vartype = "calhi_"
elif args.mcslo:
    vartype = "mcslo_"
elif args.mcshi:
    vartype = "mcshi_"
else: # Nominal
    vartype = ""

if args.loose:
    vartype += "loose_"

if args.cut_signal_box:
    vartype += "nosignal_"
    

f = args.filename

if not args.filename.endswith("evt.df"):
    print("Must run on spectrum file!")
    sys.exit(1)

spectrum_file = f.replace("evt.df", "%sevt_spectrum.df" % vartype)
if args.output_dir:
    spectrum_file = os.path.join(args.output_dir, os.path.basename(spectrum_file))
print("Will save to: %s" % spectrum_file)

if os.path.isfile(spectrum_file) and not args.overwrite:
    print("File already exists! Exiting normally.")
    sys.exit(0)

try:
    mchdf = pd.read_hdf(f, key="mch")
except:
    mchdf = pd.DataFrame()
if mchdf.index.nlevels > 2:
    mchdf.index = mchdf.index.droplevel(2)

hdr = pd.read_hdf(f, key="hdr")

issignal = mchdf.shape[0] > 0 # MeVPrtl MC

# Load the reconstructed events
evtdf = data.mc_dataset(f, "evt", mcnukey="mcnuwgt", isscalar=issignal).df.copy() if not args.data else data.onbeam_dataset(f, "evt").df.copy()

# MeVPrtl truth information
# Simulation parameters
if issignal:
    p_E = mchdf.E
    p_E.name = "p_E"

    p_M = mchdf.M
    p_M.name = "p_M"

    p_th = mchdf.C1
    p_th.name = "p_th"

    p_dist = dmagdf(mchdf.start, mchdf.enter)
    p_dist.name = "p_dist"

    p_raydist = dmagdf(mchdf.enter, mchdf.exit)
    p_raydist.name = "p_raydist"

    p_decaylength = mchdf.decay_length
    p_decaylength.name = "p_decaylength"

    evtdf[p_E.name] = p_E
    evtdf[p_M.name] = p_M
    evtdf[p_th.name] = p_th
    evtdf[p_raydist.name] = p_raydist
    evtdf[p_dist.name] = p_dist
    evtdf[p_decaylength.name] = p_decaylength
    
    evtdf["higgs"] = True
    
else: # Neutrino MC or data
    evtdf["p_E"] = np.nan
    evtdf["p_M"] = np.nan
    evtdf["p_th"] = np.nan
    evtdf["p_dist"] = np.nan
    evtdf["p_raydist"] = np.nan
    evtdf["p_dist"] = np.nan
    evtdf["p_decaylength"] = np.nan
    
    evtdf["higgs"] = False

mcdf = pd.read_hdf(f, key="mcnu")
def cohlike(df):
    nmu = df.iscc & (np.abs(df.pdg) == 14)
    npi = df.npi
    ns = df.nsm + df.nsp

    is_coh_like = (nmu + npi + ns >= 2) & (df.max_proton_ke < 0.05)

    is_coh_like = is_coh_like & InFV(df.position, 50)

    return is_coh_like

if not args.data:
    evtdf = evtdf.join(cohlike(mcdf).groupby(level=[0,1]).any().rename(("cohlike", "", "", "", "", "")))
else:
    evtdf["cohlike"] = False


if not args.data:
    pot = np.sum(hdr.pot * hdr.first_in_subrun) if not args.fix_pot else hdr.shape[0] / NU_PER_POT_BASE
    GOAL_POT = 2.41e20 
    evtdf["scale"] = evtdf.wgt.cv*GOAL_POT / pot
    evtdf["cvscale"] = evtdf.wgt.cv
else:
    pot = 0
    evtdf["scale"] = 1
    evtdf["cvscale"] = 1

if issignal and "time" in mchdf.columns: # apply timing fix for scalars
    wgt = timewgt.reweight(mchdf.time/1e3)
    evtdf["timewgt"] = wgt
else:
    evtdf["timewgt"] = 1

beam2det = np.array([ [0.921035925, 0.022715103, 0.388814672], [0, 0.998297825, -0.058321970], [-0.389477631, 0.053716629, 0.919468161]])
beamorigin = np.array([4.503730e2, 80.153901e2, 795.112945e2])
MUON_MASS = 0.1057 
PION_MASS = 0.13457
BEAMDIR = beam2det.dot(beamorigin) / np.linalg.norm(beam2det.dot(beamorigin)) 

def recop(trk):
    p = trk.rangeP.p_muon
    p[~TrkInFV(trk.end)] = trk.mcsP.fwdP_muon
    return p

scalar_mom = evtdf.trunk.trk.dir.multiply(recop(evtdf.trunk.trk), axis=0) + \
           evtdf.branch.trk.dir.multiply(recop(evtdf.branch.trk), axis=0)
    
scalar_dir = scalar_mom.divide(magdf(scalar_mom), axis=0)

beamangle = np.arccos(scalar_dir.x*BEAMDIR[0] + scalar_dir.y*BEAMDIR[1] + scalar_dir.z*BEAMDIR[2])*180/np.pi

trunk_beamcos = evtdf.trunk.trk.dir.x*BEAMDIR[0] + evtdf.trunk.trk.dir.y*BEAMDIR[1] + evtdf.trunk.trk.dir.z*BEAMDIR[2]
trunk_E = np.sqrt(recop(evtdf.trunk.trk)**2 + MUON_MASS**2)
branch_beamcos = evtdf.branch.trk.dir.x*BEAMDIR[0] + evtdf.branch.trk.dir.y*BEAMDIR[1] + evtdf.branch.trk.dir.z*BEAMDIR[2]
branch_E = np.sqrt(recop(evtdf.branch.trk)**2 + PION_MASS**2)
coh_t = (np.sin(beamangle)*magdf(scalar_mom))**2 + (trunk_E - recop(evtdf.trunk.trk)*trunk_beamcos + branch_E - recop(evtdf.branch.trk)*branch_beamcos)**2

openangle = np.arccos(dotdf(evtdf.trunk.trk.dir, evtdf.branch.trk.dir))*180/np.pi

# OPTIMIZED
MAX_SHW_LEN_CUT = 15
MAX_OTHR_TRK_LEN_CUT = 15
OTHR_CHI2U_CUT = 45
OTHR_CHI2P_CUT = 90

TRUNK_LEN_CUT = 25
BRANCH_LEN_CUT = 25
TRUNK_CHI2U_CUT = 10 # 15?
BRANCH_CHI2U_CUT = 10
CHI2P_CUT = 95
MCS_RANGE_COMP = 1.8

STUB_5mm_CUT = 20
STUB_10mm_CUT = 30
STUB_20mm_CUT = 35
STUB_30mm_CUT = 45

OPENANGLE_CUT = 70
UU_ANGLE_CUT = 5 if not args.loose else 15

# MUON ID VARIATIONS
if args.gainlo:
    chi2u_trunk = evtdf.trunk.trk.chi2pid.I2_glo.chi2u
    chi2u_branch = evtdf.branch.trk.chi2pid.I2_glo.chi2u
    chi2p_trunk = evtdf.trunk.trk.chi2pid.I2_glo.chi2p
    chi2p_branch = evtdf.branch.trk.chi2pid.I2_glo.chi2p
elif args.gainhi:
    chi2u_trunk = evtdf.trunk.trk.chi2pid.I2_ghi.chi2u
    chi2u_branch = evtdf.branch.trk.chi2pid.I2_ghi.chi2u
    chi2p_trunk = evtdf.trunk.trk.chi2pid.I2_ghi.chi2p
    chi2p_branch = evtdf.branch.trk.chi2pid.I2_ghi.chi2p
elif args.callo:
    chi2u_trunk = evtdf.trunk.trk.chi2pid.I2_callo.chi2u
    chi2u_branch = evtdf.branch.trk.chi2pid.I2_callo.chi2u
    chi2p_trunk = evtdf.trunk.trk.chi2pid.I2_callo.chi2p
    chi2p_branch = evtdf.branch.trk.chi2pid.I2_callo.chi2p
elif args.calhi:
    chi2u_trunk = evtdf.trunk.trk.chi2pid.I2_calhi.chi2u
    chi2u_branch = evtdf.branch.trk.chi2pid.I2_calhi.chi2u
    chi2p_trunk = evtdf.trunk.trk.chi2pid.I2_calhi.chi2p
    chi2p_branch = evtdf.branch.trk.chi2pid.I2_calhi.chi2p
elif not args.data: # Nominal
    chi2u_trunk = evtdf.trunk.trk.chi2pid.I2.chi2u
    chi2u_branch = evtdf.branch.trk.chi2pid.I2.chi2u
    chi2p_trunk = evtdf.trunk.trk.chi2pid.I2.chi2p
    chi2p_branch = evtdf.branch.trk.chi2pid.I2.chi2p
else: # different VAR name for data
    chi2u_trunk = evtdf.trunk.trk.chi2pid.I2.chi2_muon
    chi2u_branch = evtdf.branch.trk.chi2pid.I2.chi2_muon
    chi2p_trunk = evtdf.trunk.trk.chi2pid.I2.chi2_proton
    chi2p_branch = evtdf.branch.trk.chi2pid.I2.chi2_proton

trunk_mcs_range = evtdf.trunk.trk.mcsP.fwdP_muon / evtdf.trunk.trk.rangeP.p_muon
branch_mcs_range = evtdf.branch.trk.mcsP.fwdP_muon / evtdf.branch.trk.rangeP.p_muon

if args.mcslo:
    trunk_mcs_range -= 0.03
    branch_mcs_range -= 0.03
elif args.mcshi:
    trunk_mcs_range += 0.03
    branch_mcs_range += 0.03
else: # Nominal
    pass

# STUB ID VARIATIONS
if args.callo and not issignal: # only apply stub variations to neutrinos
    stubl0_5_dedx = evtdf.stub.l0_5cm.dedx_callo*evtdf.stub.l0_5cm.length
    stubl1_dedx = evtdf.stub.l1cm.dedx_callo*evtdf.stub.l1cm.length
    stubl2_dedx = evtdf.stub.l2cm.dedx_callo*evtdf.stub.l2cm.length
    stubl3_dedx = evtdf.stub.l3cm.dedx_callo*evtdf.stub.l3cm.length
elif args.calhi and not issignal:
    stubl0_5_dedx = evtdf.stub.l0_5cm.dedx_calhi*evtdf.stub.l0_5cm.length
    stubl1_dedx = evtdf.stub.l1cm.dedx_calhi*evtdf.stub.l1cm.length
    stubl2_dedx = evtdf.stub.l2cm.dedx_calhi*evtdf.stub.l2cm.length
    stubl3_dedx = evtdf.stub.l3cm.dedx_calhi*evtdf.stub.l3cm.length
else: # Nominal
    stubl0_5_dedx = evtdf.stub.l0_5cm.dedx*evtdf.stub.l0_5cm.length
    stubl1_dedx = evtdf.stub.l1cm.dedx*evtdf.stub.l1cm.length
    stubl2_dedx = evtdf.stub.l2cm.dedx*evtdf.stub.l2cm.length
    stubl3_dedx = evtdf.stub.l3cm.dedx*evtdf.stub.l3cm.length

preselect = [
    SlcInFV(evtdf.slc.vertex),
    TrkInFV(evtdf.trunk.trk.end),
    TrkInFV(evtdf.branch.trk.end)
]

objs = [
    np.isnan(evtdf.max_shw_len_primary) | (evtdf.max_shw_len_primary < MAX_SHW_LEN_CUT),
    np.isnan(evtdf.max_othr_trk_len) | (evtdf.max_othr_trk_len < MAX_OTHR_TRK_LEN_CUT),
    np.isnan(evtdf.max_othr_chi2_muon) | (evtdf.max_othr_chi2_muon < OTHR_CHI2U_CUT),
    np.isnan(evtdf.min_othr_chi2_proton_all) | (evtdf.min_othr_chi2_proton_all > OTHR_CHI2P_CUT),
    
]

hasstub = \
    (stubl0_5_dedx > STUB_5mm_CUT) |\
    (stubl1_dedx > STUB_10mm_CUT) |\
    (stubl2_dedx > STUB_20mm_CUT) |\
    (stubl3_dedx > STUB_30mm_CUT)
    
muid = [
    evtdf.trunk.trk.len > TRUNK_LEN_CUT,
    evtdf.branch.trk.len > BRANCH_LEN_CUT,
    chi2u_trunk < TRUNK_CHI2U_CUT,
    chi2u_branch < BRANCH_CHI2U_CUT, 
    chi2p_trunk > CHI2P_CUT,
    chi2p_branch > CHI2P_CUT,
    trunk_mcs_range < MCS_RANGE_COMP,
    branch_mcs_range < MCS_RANGE_COMP,
]

kinematics = [
    openangle < OPENANGLE_CUT,
    beamangle < UU_ANGLE_CUT,
]

if args.cut_signal_box:
    kinematics[-1] = kinematics[-1] & (beamangle > 5)

stubs = [
    ~hasstub,
]

def cut_all(l):
    l0 = l[0]
    for c in l[1:]:
        l0 = l0 & c
    return l0

cuts = preselect + objs + muid + stubs + kinematics

passcut = cut_all(cuts)

passcut_nokinematics = cut_all(preselect + muid + stubs)

Etrunk = np.sqrt(MUON_MASS**2 + recop(evtdf.trunk.trk)**2)
Ebranch = np.sqrt(MUON_MASS**2 + recop(evtdf.branch.trk)**2)

s_mass = np.sqrt(2)*np.sqrt(MUON_MASS**2 + Etrunk*Ebranch - recop(evtdf.trunk.trk)*recop(evtdf.branch.trk)*dotdf(evtdf.branch.trk.dir, evtdf.trunk.trk.dir))

# Compute track splitting variables and apply corrections
if args.correct_split_tracks and not args.data:
     import track_splitting
     cathode_cross_weight, cathode_cross_weight_err, gap_cross_weight, gap_cross_weight_err = track_splitting.correct(evtdf)
else:
    cathode_cross_weight = pd.Series(1, evtdf.index)
    cathode_cross_weight_err = pd.Series(0, evtdf.index)
    gap_cross_weight = pd.Series(1, evtdf.index)
    gap_cross_weight_err = pd.Series(0, evtdf.index)

# Update weights
evtdf.scale = evtdf.scale*cathode_cross_weight*gap_cross_weight
evtdf.cvscale = evtdf.cvscale*cathode_cross_weight*gap_cross_weight

spectradf = pd.DataFrame(
    data={
        "mass": s_mass[passcut],
        "openangle": openangle[passcut],
        "thnumi": beamangle[passcut], 
        "tuu": coh_t[passcut],
        "chi2u_trunk": chi2u_trunk[passcut], 
        "chi2u_branch": chi2u_branch[passcut], 
        "chi2p_trunk": chi2p_trunk[passcut], 
        "chi2p_branch": chi2p_branch[passcut], 
        "signal": evtdf.higgs[passcut],
        "cohlike": evtdf.cohlike[passcut],
        "weight": evtdf.scale[passcut],
        "pot": pot,
        "cvweight": evtdf.cvscale[passcut],
        "timeweight": evtdf.timewgt[passcut],
        "cathode_weight": cathode_cross_weight[passcut],
        "cathode_weight_err": cathode_cross_weight_err[passcut],
        "gap_weight": gap_cross_weight[passcut],
        "gap_weight_err": gap_cross_weight_err[passcut],
        "genie_mode": evtdf.slc.truth.genie_mode[passcut],
        "nu_iscc": evtdf.slc.truth.iscc[passcut],
        "nu_pdg": evtdf.slc.truth.pdg[passcut],
        "nu_E": evtdf.slc.truth.E[passcut],
        "nu_Q2": evtdf.slc.truth.Q2[passcut],
        "nu_W": evtdf.slc.truth.w[passcut],
        "p_E": evtdf.p_E[passcut],
        "p_M": evtdf.p_M[passcut],
        "p_th": evtdf.p_th[passcut],
        "p_dist": evtdf.p_dist[passcut],
        "p_raydist": evtdf.p_raydist[passcut],
        "p_decaylength": evtdf.p_decaylength[passcut],
    },
)

# systematic uncertainties
if not args.data:
    systdfs = []
    for syst in ["g4", "xsec", "flux", "coh", "all"]:
        systdf = evtdf.wgt[syst][passcut]
        systdf.columns = [("wgt_%s_" % syst) + "_".join([r for r in c if r]) for c in systdf.columns]
        systdfs.append(systdf)

    systdf = pd.concat(systdfs, axis=1)

    # Add in weights
    spectradf = pd.concat([spectradf, systdf], axis=1).reset_index(drop=True)

if args.remove_cohlike:
    spectradf = spectradf[~spectradf.cohlike]

# Save
spectradf.to_hdf(spectrum_file, key="s")
