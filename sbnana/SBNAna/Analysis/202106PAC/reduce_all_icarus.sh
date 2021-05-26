#!/bin/bash

# This is the script used to apply reduce_icarus.C to all the Icarus CAFs. This
# takes a very long time, so we (ab)use a 64 core machine. This is definitely
# not the best way to write this macro, but it works well enough.

if [ `hostname -s` != sbndbuild02 ]
then
    echo This script is explicitly designed to run \(cheekily\) on sbndbuild02. It runs 64 processes in parallel. If you run it on a machine with fewer than 64 cores you are likely to have a bad time.
    exit 1
fi

INDIR=/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/
OUTDIR=/pnfs/sbn/persistent/analysis/CAF/202106PAC/workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus/reduced_cafs/

MACRODIR=`pwd`

for vi in '0 1 2 3' '4 5 6 7' '8 9 a b' 'c d e f'
do
    for i in $vi
    do
        for j in `seq 0 9` a b c d e f
        do
            # ROOT generates temporary files in the same directory as the
            # macro, which makes a big mess and fails when 64 copies do it all
            # at once. Have each process work in a seperate directory (and
            # never clean them up... naughty).
            cd `mktemp -d`
            cp $MACRODIR/reduce_icarus.C .
            cafe -bq reduce_icarus.C $INDIR/$i/$j/'*.root' $OUTDIR/${i}_${j}_x.caf.root &
        done
    done
    wait # Wait for all 64 jobs to finish before starting the next batch
done
