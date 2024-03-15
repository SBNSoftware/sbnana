import argparse
import pandas as pd
import numpy as np
import sys

# configured
DETSYST_NU = 15.426
DETSYST_DCY = 9.8711
NUMISYST = 8.6

NUNI = 100
np.random.seed(103) # It's October 3rd
mass_bins = np.array([0.21, 0.225, 0.245, 0.265, 0.29, 0.32, 0.35, 0.4, 0.5])

parser = argparse.ArgumentParser(
    prog="add_detsyst_2spectrum.py",
    description="Add detector uncertainties to spectrum dataframe.",
)
parser.add_argument("filename")

args = parser.parse_args()
f = args.filename

if not f.endswith("spectrum.df"):
    print("Must run on spectrum file!")
    sys.exit(1)

spectrum_file = f.replace("spectrum.df", "spectrum_detsyst.df")
print("Will save to: %s" % spectrum_file)

DETSYST = DETSYST_DCY if "Higgs" in spectrum_file else DETSYST_NU

df = pd.read_hdf(f, key="s")

# Apply time-weighting
if "timeweight" in df.columns:
    df.weight *= df.timeweight
    df.cvweight *= df.timeweight

# Detector systematic uncertainty (normalization)
detsyst_univ = 1 + np.random.normal(size=NUNI)*DETSYST/100

numisyst_univ = 1 + np.random.normal(size=NUNI)*NUMISYST/100

# Cathode crossing uncertainty
cathode_univ = np.random.normal(size=NUNI)

# Induction gap uncertainty
gap_univ = np.random.normal(size=NUNI)

# M.C. stat uncertainty
mcstat_uncs = []
for mlo, mhi in zip(mass_bins[:-1], mass_bins[1:]):
    when = (df.mass > mlo) & (df.mass < mhi) & (df.thnumi < 5)
    mcstat_uncs.append(np.nan_to_num(np.sqrt((df.weight[when]**2).sum())/df.weight[when].sum()))
mcstat = pd.cut(df.mass, mass_bins, labels=mcstat_uncs, ordered=False).astype(float).fillna(0)

mcstat_univ = np.random.normal(size=(NUNI, len(mass_bins)-1)) 

for i in range(NUNI):
    # multiply the normalization, cathode-cross, and induction-gap-cross uncertainties
    cathode_err_wgt = np.maximum(1 + cathode_univ[i]*df.cathode_weight_err, 0)
    gap_err_wgt = np.maximum(1 + gap_univ[i]*df.gap_weight_err*df.gap_weight_err, 0)
    df["wgt_det_univ_%i" % i] = np.maximum(detsyst_univ[i]*cathode_err_wgt*gap_err_wgt, 0)
    df["wgt_all_univ_%i" % i] *= df["wgt_det_univ_%i" % i]

    mcstat_univ_var = pd.cut(df.mass, mass_bins, labels=mcstat_univ[i], ordered=False).astype(float).fillna(0)
    mcstat_univ_var[df.thnumi > 5] = 0 # assume ~0 MC stat uncertainty outside signal box 
    df["wgt_mcstat_univ_%i" % i] = np.maximum(1 + mcstat_univ_var * mcstat, 0)
    df["wgt_all_univ_%i" % i] *= df["wgt_mcstat_univ_%i" % i]

    # additional flux uncertainty
    df["wgt_flux_univ_%i" % i] *= numisyst_univ[i]
    df["wgt_all_univ_%i" % i] *= numisyst_univ[i]

    df = df.copy() # undo fragmentation

# Save
df.to_hdf(spectrum_file, key="s")
