source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup python v3_9_2
setup xrootd v5_1_0 -q e20:p392:prof
setup hdf5 v1_12_0a -q e20:prof
setup sbnana v09_37_02_01 -q e20:prof

python -m venv env
source env/bin/activate
which python
pip install --upgrade pip
pip install wheel setuptools
pip install -r requirements.txt

# Needed to install xrootd -- which, by the way, is super annoying
OLDPATH=$PATH
PATH=$PATH:$PWD
ln -s /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_22_2/Linux64bit+3.10-2.17/bin/cmake cmake3
which cmake3
wget https://files.pythonhosted.org/packages/2c/91/c03d649236b3af7720b2ea5561abfbbf3baed3722da485cfea07e565cf57/xrootd-5.5.4.tar.gz
tar -zxvf xrootd-5.5.4.tar.gz 
rm xrootd-5.5.4.tar.gz
cd xrootd-5.5.4/
python setup.py install
cd ..
PATH=$OLDPATH
