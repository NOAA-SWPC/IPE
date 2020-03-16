#!/bin/bash

## USER SETUP
# IPE_State.apex.*.h5 location
IPE_INPUT_DIRECTORY=/gpfs/dell2/ptmp/George.Millward/IPEstandalone_weimer_test/
# where you want the IPE_Params.geo.*.nc4 files
IPE_OUTPUT_DIRECTORY=${IPE_INPUT_DIRECTORY}/netcdf
###### USUALLY NO NEED TO CHANGE ANYTHING BELOW THIS LINE ######
if [[ -e /scratch1 ]] ; then
  ## HERA GRID LOCATION
  IPE_GRID=/scratch1/NCEPDEV/swpc/WAM-IPE_DATA/IPE_FIX/IPE_Grid.h5
  ## MODULEFILE LOAD
  module use -a /contrib/modulefiles
  module load anaconda
elif [[ -e /gpfs/dell2 ]] ; then
  ## WCOSS DELL P3 GRID LOCATION
  IPE_GRID=/gpfs/dell2/swpc/noscrub/George.Millward/ipe_grid/IPE_Grid.h5
  ## MODULEFILE LOAD

#  module load python/2.7.14
  module load python/2.7.13
fi
## DIRECTORY SETUP
mkdir -p $IPE_OUTPUT_DIRECTORY
rm ${IPE_OUTPUT_DIRECTORY}/*.nc4
## RUN
python /gpfs/dell2/swpc/noscrub/George.Millward/gsm_repo2/GSMWAM-IPE/IPE/scripts/post_minimal.py -g $IPE_GRID -i $IPE_INPUT_DIRECTORY -o $IPE_OUTPUT_DIRECTORY
