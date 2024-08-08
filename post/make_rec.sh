#!/bin/bash

module purge
module load releases/2022b
module load ELIC_Python/1-foss-2022b

#-------------
# Parameters |
#-------------
exp_name='2018-2025_all_vars_ANN_UKESM'
outfolder_rec='/cyfast/dalaiden/20th_reconstruction_hgs_fogt/DA_exps_outputs'
#-------------

echo '-----------------------'
echo "Reconstruct the fields|"
echo '-----------------------'

export PYTHONWARNINGS="ignore"
python -W ignore make_prior.py $exp_name
python -W ignore make_posterior.py $exp_name $outfolder_rec
