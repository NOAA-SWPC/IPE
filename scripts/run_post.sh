#!/bin/bash

## USER SETUP
# IPE_State.apex.*.h5 location
IPE_INPUT_DIRECTORY=/scratch4/NCEPDEV/stmp3/Adam.Kubaryk/comp_hdf5_grid
# where you want the IPE_Params.geo.*.nc4 files
IPE_OUTPUT_DIRECTORY=$IPE_INPUT_DIRECTORY/netcdf

###### USUALLY NO NEED TO CHANGE ANYTHING BELOW THIS LINE ######

## THEIA GRID LOCATION
IPE_GRID=/scratch3/NCEPDEV/swpc/noscrub/refactored_ipe_input_decks/IPE_Grid.h5
## MODULEFILE LOAD
module purge
module use -a ../../NEMS/src/conf
module use -a /contrib/modulefiles
module load modules.nems
module load anaconda
## DIRECTORY SETUP
mkdir -p $IPE_OUTPUT_DIRECTORY
## RUN
python post.py -g $IPE_GRID -i $IPE_INPUT_DIRECTORY -o $IPE_OUTPUT_DIRECTORY
