#source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup python v3_9_2
setup xrootd v5_1_0 -q e20:p392:prof
setup hdf5 v1_12_0a -q e20:prof
setup sbnana v09_69_01 -q e20:prof

# Needed to install xrootd
OLDPATH=$PATH
PATH=$PATH:$PWD
ln -s /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_22_2/Linux64bit+3.10-2.17/bin/cmake cmake3
which cmake3

python -m venv env
source env/bin/activate
which python
pip install --upgrade pip
pip install wheel setuptools
pip install -r requirements.txt

PATH=$OLDPATH
