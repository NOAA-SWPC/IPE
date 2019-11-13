#!/bin/bash

## USER SETUP
# IPE_State.apex.*.h5 location
IPE_INPUT_DIRECTORY=/scratch1/NCEPDEV/stmp2/Tzu-Wei.Fang/GSMWAM-IPE_dynamo
# where you want the IPE_Params.geo.*.nc4 files
#IPE_OUTPUT_DIRECTORY=$IPE_INPUT_DIRECTORY/netcdf
IPE_OUTPUT_DIRECTORY=/scratch1/NCEPDEV/stmp4/Adam.Kubaryk/tec_geo_test
###### USUALLY NO NEED TO CHANGE ANYTHING BELOW THIS LINE ######

## HERA GRID LOCATION
IPE_GRID=/scratch1/NCEPDEV/swpc/WAM-IPE_DATA/IPE_FIX/IPE_Grid.h5
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
