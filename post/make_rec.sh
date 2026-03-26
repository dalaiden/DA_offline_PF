#!/bin/bash

module purge
module load releases/2025b
module load ELIC_Python/1-foss-2025b

#-------------
# Parameters |
#-------------
exp_name='1800-2020_all_records_CESM1-LM_20240725_1.5std_v2'
outfolder_rec='/cyfast/hxue/DA_output'
#-------------

echo '-----------------------'
echo "Reconstruct the fields|"
echo '-----------------------'

export PYTHONWARNINGS="ignore"
python -W ignore make_prior_v2.py $exp_name
python -W ignore make_posterior_v2.py $exp_name $outfolder_rec
