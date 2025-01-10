OUTDIR=/icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/

python make_spectrum_df.py -od $OUTDIR -cst -rc /icarus/data/users/gputnam/DMCP2023G/mc-F/F-MCNuPhase2_evt.df
python make_spectrum_df.py -od $OUTDIR --loose -cst -rc /icarus/data/users/gputnam/DMCP2023G/mc-F/F-MCNuPhase2_evt.df

for arg in gl gh cl ch ml mh
do
  python make_spectrum_df.py -od $OUTDIR --loose -cst -${arg} /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_nom_evt.df &
  python make_spectrum_df.py -od $OUTDIR -cst -${arg} /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike_nom_evt.df &
  wait
done

for f in /icarus/data/users/gputnam/DMCP2023G/mc-F/F-CohLike*evt.df
do
  python make_spectrum_df.py -od $OUTDIR -cst --loose $f &
  python make_spectrum_df.py -od $OUTDIR -cst $f &
  wait
done

for f in `ls /icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom*evt.df | grep -v polaris | grep -v nowgt`
do
  python make_spectrum_df.py -od $OUTDIR -cst --loose $f &
  python make_spectrum_df.py -od $OUTDIR -cst $f &
  wait
done

for f in `ls /icarus/data/users/gputnam/DMCP2023G/mc-FB/FB-CohLike_nom*evt.df | grep polaris | grep -v nowgt`
do
  python make_spectrum_df.py -fp -od $OUTDIR -cst --loose $f &
  python make_spectrum_df.py -fp -od $OUTDIR -cst $f &
done
wait

python add_dfs_fixpot.py /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-CohLike_nom-all_evt_spectrum.df   /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-CohLike_nom_evt_spectrum.df `ls /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/FB-CohLike_nom-* | grep -v loose`
python add_dfs_fixpot.py /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-CohLike_nom-all_loose_evt_spectrum.df   /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-CohLike_nom_loose_evt_spectrum.df `ls /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/FB-CohLike_nom-* | grep loose`

python add_dfs.py /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-Nu_loose_evt_spectrum.df /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-CohLike_nom-all_loose_evt_spectrum.df /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-MCNuPhase2_loose_evt_spectrum.df 
python add_dfs.py /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-Nu_evt_spectrum.df /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-CohLike_nom-all_evt_spectrum.df /icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/F-MCNuPhase2_evt_spectrum.df 

for M in 220 240 260 280 300 320 340
do
  for f in `ls /icarus/data/users/gputnam/DMCP2023G/mc-F/F2-Higgs_M*_evt.df | grep $M | grep -v _nom_`
  do
    python make_spectrum_df.py -od $OUTDIR -cst $f &
  done
  wait
done

for f in  /icarus/data/users/gputnam/DMCP2023G/mc-F-iter2/F2-Higgs_M*_nom_evt.df
do
  python make_spectrum_df.py -od $OUTDIR -cst $f &
  python make_spectrum_df.py -od $OUTDIR -cst --loose $f &
  for arg in gl gh cl ch ml mh
  do
    python make_spectrum_df.py -od $OUTDIR -cst -${arg} $f &
  done
  wait
done
