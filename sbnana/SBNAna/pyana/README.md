## To initialize environment

source init.sh

## To setup environment in new terminal after initialization (only need to do this once after initialization)


source setup.sh

## Generate a HDF5 file from a .flat.caf.root file

`python run.py -c configs/config.py -o /path/to/output/data.df -i /path/to/inputs/\*.flat.caf.root`

## Generate a HDF5 file from a .list file
.list contains line by line list of flat cafs

`python run.py -c configs/config.py -o /path/to/output/data.df -i /path/to/inputs/\*.list`

## Run jupyter notebooks in nb/

jupyter notebook --no-browser
