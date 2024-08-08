# Offline Data Assimilation (particle filter; UCLouvain-ELIC)

This is the offline data assimilation code utilizing a particle filter, as outlined in Dubinkina et al. (2011).  The program has been employed in several recent publications. For further details about the method, please refer to these publications. The repository contains the data assimilation framework developped for reconstructing the Antarctic sea ice and related variables over 1958-2023 using station-based observations as described in Goosse et al. (2024; ). For further uses, please contact Quentin Dalaiden ([quentin.dalaiden@uclouvain.be](quentin.dalaiden@uclouvain.be)) to adapt the framework. 

## Installation

- git clone https://github.com/dalaiden/DA_offline_PF.git

From the folder main folder:
- Open the file `src/modules.load` and check that you have all requested libraries at the right place (currently templates are available for the clusters `storm`, `coriolis` and `cyclone`).
- Execute the Makefile: `cd src; . modules.load; make`.
- Check that you have your executables: `src/moddata_co` and `src/PartFilter`.

## Program

`./assim_offline.sc` is the main script calling the data assimilation. Information about the main parameters as well as the progress of the program will appear in the terminal. At the end of the assimilation, output files containing the weights of each ensemble member (i.e., particle) are created in the folder `rundir/<exp_name>/output_fcost`. Once finished, the posterior is created and exported to a defined folder set by the user.

## Experiment preparation

### Input files

The program requires several input files: 

- the observations to assimilate
- the model ensemble member(s) to build the prior
- the file containing the error associated with the assimilated observations

The python script `make_inputs.py` creates all these required files.

Once the script finished, the path of the files must be specified in `moddata_co.namelist_brut`, except the model ensemble members, which are included in `assim_offline.sc`:

**Line 49-76** Template showing the order in which the information about the variables assimilated should be defined between the lines 65 to 80. One bloc per variable assimilated should be filled in, with always the same structure:

- var name (mandatory). The name of the variable to be assimilated, as specified in the netcdf files.
- var type (mandatory). The realm of the variable assimilated, ie. atmos or atmos{2..6}. When several atmopsheric variables are assimilated, add a number after *atmos* (e.g., *atmos2*)
- 4D structure (mandatory). The dimension of the netcdf file. Do not change.
- weight (mandatory). The weight that is given, when more than one variable is assimilated, to each variable assimilated. Do not change.
- sqrt(Ci) (mandatory). Not used anymore.
- Sigma (mandatory). The weight for the model covariance matrix. 0 if cov matrix should be diagonal. 
- Number of domains (mandatory) Number of domains used. Set to 1 if no spatial averages are performed.
- Box Type (optional). Not used anymore
- ObsFileStart: startDateY (mandatory). Must be 1. Do not change.
- ObsFileStart: startDateD (mandatory). Must be 1. Do not change.
- DataObs (mandatory). Netcdf file including the data that should be assimilated, with full path.
- RefObs (mandatory). Netcdf file of the reference file of the data assimilated, with full path.
- ErrObs (mandatory). Netcdf file containing the error associated with the data assimilated, with full path.
- RefMod (mandatory). The reference file of the model simulations, with full path.
- Datacov (optional) Not used anymore.

### Experimental design

The timing of the experiment and the options related to the ensemble size are specified in the file assim_offline.sc. For normal use, you are not supposed to change anything beyond the line 51. Here are the parameters to edit:

- **exp_name** is the name of your experiment. It cannot contain any space or special character. A folder of that name will be created in the directory rundir, containing the output of the experiment.
- **duration_model** is the length of the model simulations, expressed in months. It must be a multiple of frequence_assim.
- **duration_data** is the length of the data file that is going to be assimilated. It can be different than duration_model. If duration_data > duration_model, increase_ensemble_size has to be equal to "1".
- **frequence_assim** is the frequency of the assimilation. It can be equal to 1 (monthly assimilation), 12 (annual mean assimilation), or higher. In the latter case, it has to be a multiple of 12. If frequence_assim > 12, it is not yet possible to increase the ensemble size (the option increase_ensemble_size will have to be equal to 0.)
- **increase_ensemble_size** determines whether the filter artificially increases the ensemble size by selecting additional particles with a date different from the one of the observations. To activate this option, set increase_ensemble_size to 1. To respect the right timing by only taking into account the particles of the year corresponding to the one in the data file, set increase_ensemble_size to 0.
- **frequence_sampling** is the frequency of the sampling when you ask the filter to select other particles to the ones of the data (increase_ensemble_size="1"). If the frequency of the assimilation is annual frequence_assim="12", you can set frequence_sampling to whatever value you like. For instance, if frequence_sampling="5", the filter will also consider, besides the particles of the actual year, particles of 1 year out of 5 over the whole period of the assimilation. If the assimilation is monthly (frequence_assim="1"), frequence_sampling should be a multiple of 12, in order to select only particles of the appropriate month. If increase_ensemble_size="0", this parameter won’t be used.
- **directory_input** The last variable to edit in this file is directory_input. It is the full path of the folder containing the model ensemble members, without the last /. This folder should not contain any other netcdf files. The name of the different nc files
should only differ from the numbers of the ensemble members. For instance, the files can be named: file_1.nc, file_2.nc, etc. directory_input also contains the realm of the variable to be assimilated, ie. atmos for variables belonging to the atmospheric grid and ocean for variables belonging to the oceanic one., separated by \<space>:\<space>. 
For instance:
`declare -a directory_input=("/address_atmos_simulations : atmos")`
When two (or more) variables are assimilated, here is how it should be written:
`declare -a directory_input=("/address_atmos_simulations : atmos" "/address_ocean_simulations : ocean")`

## Launch an experiment and get results

Execute assim_offline.sc. Information about the progress of the program appears on the screen.

Once the first assimilation is finished, an output file containing the the fcosts of the different particles is created in the folder `rundir/<exp_name>/output_fcost`. There is one fcost file per assimilation, with one line per particle.

If the program crashes or if you want more information, you can have a look at the log files, located in the folder `rundir/<exp_name>/output_log`. There is one log file per particle.