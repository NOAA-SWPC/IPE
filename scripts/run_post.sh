#!/bin/bash

## USER SETUP
# IPE_State.apex.*.h5 location
IPE_INPUT_DIRECTORY=/scratch1/NCEPDEV/stmp2/Tzu-Wei.Fang/GSMWAM-IPE_dynamo
# where you want the IPE_Params.geo.*.nc4 files
IPE_OUTPUT_DIRECTORY=$IPE_INPUT_DIRECTORY/netcdf
###### USUALLY NO NEED TO CHANGE ANYTHING BELOW THIS LINE ######
if [[ -e /scratch1 ]] ; then
  ## HERA GRID LOCATION
  IPE_GRID=/scratch1/NCEPDEV/swpc/WAM-IPE_DATA/IPE_FIX/IPE_Grid.h5
  ## MODULEFILE LOAD
  module use -a /contrib/modulefiles
  module load anaconda
elif [[ -e /gpfs/dell2 ]] ; then
  ## WCOSS DELL P3 GRID LOCATION
  IPE_GRID=/gpfs/dell2/emc/modeling/noscrub/Adam.Kubaryk/refactored_ipe_input_decks/IPE_Grid.h5
  ## MODULEFILE LOAD
  module load python/2.7.14
fi
## DIRECTORY SETUP
mkdir -p $IPE_OUTPUT_DIRECTORY
## RUN
python post.py -g $IPE_GRID -i $IPE_INPUT_DIRECTORY -o $IPE_OUTPUT_DIRECTORY
