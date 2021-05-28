#!/bin/bash

# Run this script after reduce_all_icarus.sh is complete
#
# Flattening is about 5sec per file, for 20mins total

TOPDIR=/pnfs/sbn/persistent/analysis/CAF/202106PAC/workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus/

INDIR=$TOPDIR/reduced_cafs/

OUTDIR=$TOPDIR/reduced_flatcafs/

for i in `seq 0 9` a b c d e f
do
    for j in `seq 0 9` a b c d e f
    do
        flatten_caf $INDIR/${i}_${j}_x.caf.root $OUTDIR/${i}_${j}_x.flat.caf.root
        # Protect against accidents
        chmod -w $OUTDIR/${i}_${j}_x.flat.caf.root
        # Actually, let's do the input too. The parallel way they were made
        # makes it harder to do at that step.
        chmod -w $INDIR/${i}_${j}_x.caf.root
    done
done
