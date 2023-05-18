## To initialize environment

source init.sh

## To setup environment in new terminal after initialization (only need to do this once after initialization)


source setup.sh

# Generate a HDF5 file from a .flat.caf.root file

python run.py configs/.py /path/to/output/data.df /path/to/inputs/\*.flat.caf.root

# Run jupyter notebooks in nb/

jupyter notebook --no-browser
