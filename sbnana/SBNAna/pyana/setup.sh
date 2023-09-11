#source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup python v3_9_2
setup xrootd v5_1_0 -q e20:p392:prof
setup hdf5 v1_12_0a -q e20:prof
setup sbnana v09_69_01 -q e20:prof

source env/bin/activate

export PYTHONPATH=$PYTHONPATH:$PWD
