#!/bin/bash

module purge
module load releases/2022b
module load ELIC_Python/1-foss-2022b

#-------------
# Parameters |
#-------------
exp_name='1998-2023_all_vars_1dot5std_SH_500km-grid_sx200_xy100_tassim_5'
outfolder_rec='/nas07/dalaiden/cyfast/paleoPF_ant/DA_exps_outputs'
#-------------

echo '-----------------------'
echo "Reconstruct the fields|"
echo '-----------------------'

export PYTHONWARNINGS="ignore"
python -W ignore make_prior.py $exp_name
python -W ignore make_posterior.py $exp_name $outfolder_rec
