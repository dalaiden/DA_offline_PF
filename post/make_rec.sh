#!/bin/bash

module purge
module load releases/2022b
module load ELIC_Python/1-foss-2022b

#-------------
# Parameters |
#-------------
exp_name='1850-2020_all_records_CESM1-LM_20240422_v2'
outfolder_rec='/cyfast/hxue/DA_output'
#-------------

echo '-----------------------'
echo "Reconstruct the fields|"
echo '-----------------------'

export PYTHONWARNINGS="ignore"
python -W ignore make_prior_VPD.py $exp_name
python -W ignore make_posterior_VPD.py $exp_name $outfolder_rec
