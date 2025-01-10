export machine=${HOSTNAME}
if [[ $machine == *sbnd* ]]; then
  echo "working on a sbnd machine"
  source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
fi
if [[ $machine == *icarus* ]]; then
  echo "working on a icarus machine"
  source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
fi
setup python v3_9_2
setup hdf5 v1_12_0a -q e20:prof

source env/bin/activate

export PYTHONPATH=$PYTHONPATH:$PWD/..
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VIRTUAL_ENV/lib/python3.9/site-packages/pyxrootd/lib64
