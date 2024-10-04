source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup python v3_9_2
setup hdf5 v1_12_0a -q e20:prof
setup sbnana v09_37_02_01 -q e20:prof
unsetup xrootd

source env2/bin/activate

export PYTHONPATH=$PYTHONPATH:$PWD/..
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VIRTUAL_ENV/lib/python3.9/site-packages/pyxrootd/lib64
