#source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup python v3_9_13
setup xrootd v5_4_3b -q e20:p3913:prof
setup hdf5 v1_12_0a -q e20:prof
setup sbnana v09_75_03 -q e20:prof

source env/bin/activate

export PYTHONPATH=$PYTHONPATH:$PWD
