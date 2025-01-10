source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup python v3_9_2
setup hdf5 v1_12_0a -q e20:prof
unsetup xrootd

python -m venv env	
source env/bin/activate
which python
pip install wheel setuptools
pip install -r requirements.txt

# Needed to install xrootd -- which, by the way, is super annoying
OLDPATH=$PATH
PATH=$PATH:$PWD
ln -s /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_22_2/Linux64bit+3.10-2.17/bin/cmake cmake3
which cmake3
wget https://files.pythonhosted.org/packages/e6/fd/bb238713eaede197919e12c15264f597e16302301294088af028bdff991c/xrootd-5.5.1.tar.gz
tar -zxvf xrootd-5.5.1.tar.gz 
rm xrootd-5.5.1.tar.gz
cd xrootd-5.5.1/
python setup.py install
cd ..
PATH=$OLDPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VIRTUAL_ENV/lib/python3.9/site-packages/pyxrootd/lib64
