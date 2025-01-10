PROCESSDIR=/icarus/data/users/gputnam/DMCP2023G/mc-F-spectra
python add_detsyst_2spectrum.py $PROCESSDIR/F-Nu_loose_evt_spectrum.df
for f in $PROCESSDIR/F*-Higgs_M*_loose_evt_spectrum.df
do
  python add_detsyst_2spectrum.py $f
done
