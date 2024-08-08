# Offline Data Assimilation (particle filter; UCLouvain-ELIC)

This is the offline data assimilation code utilizing a particle filter, as outlined in Dubinkina et al. (2011).  The program has been employed in several recent publications (see *References* section). For further details about the method, please refer to these publications. The repository contains the data assimilation framework developped for reconstructing the Antarctic sea ice and related variables over 1958-2023 using station-based observations as described in Goosse et al. (2024; ). For further uses, please contact Quentin Dalaiden ([quentin.dalaiden@uclouvain.be](quentin.dalaiden@uclouvain.be)) to adapt the framework. 

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
- Box Type (optional). Not used anymore.
- ObsFileStart: startDateY (mandatory). Must be 1. Do not change.
- ObsFileStart: startDateD (mandatory). Must be 1. Do not change.
- DataObs (mandatory). Netcdf file including the data that should be assimilated, with full path.
- RefObs (mandatory). Netcdf file of the reference file of the data assimilated, with full path.
- ErrObs (mandatory). Netcdf file containing the error associated with the data assimilated, with full path.
- RefMod (mandatory). The reference file of the model simulations, with full path.
- Datacov (optional) Not used anymore.

### Experimental design

The timing of the experiment and the options related to the ensemble size are specified in the file assim_offline.sc. For normal use, you are not supposed to change anything beyond the line 51. Here are the parameters to edit:

- **first_year**. The start of the assimilation.
- **exp_name**. The name of the experiment. It cannot contain any space or special character. A folder of that name will be created in the directory `rundir`, containing the output of the experiment.
- **moddata_co_file**. The file containing the paramteres of the assimilation.
- **make_posterior**. If `"True"`, the posterior is calculated at the end of the assimilation. The posterior can be produced later using the `post/make_rec.sh` script. All the variables to reconstruct are specified in `post/make_prior.py` and `post/make_posterior.py`.
- **outfolder_rec**. The directory where the posterior will be exported.
- **frequence_sampling**. The frequency of the sampling. For instance, if `frequence_sampling="1"`, the filter will consider all years of the ensemble members. If `frequence_sampling="10"`, the size of the prior will be 10 times smaller.
- **directory_input** The last variable to edit in this file is directory_input. It is the full path of the folder containing the model ensemble members, without the last */*. This folder should not contain any other netcdf files. The name of the different nc files should only differ from the numbers of the ensemble members. For instance, the files can be named: file_001.nc, file_002.nc, etc. directory_input also contains the realm of the variable to be assimilated, i.e., *atmos* or *atmos2*, separated by \<space>:\<space>. 
For instance:
`declare -a directory_input=("/address_atmos_simulations : atmos")`
When two (or more) variables are assimilated, here is how it should be written:
`declare -a directory_input=("/address_atmos_simulations : atmos" "/address_ocean_simulations : ocean")`

## Launch an experiment and get results

Execute assim_offline.sc. Information about the progress of the program appears on the screen.

Once the first assimilation is finished, an output file containing the the fcosts of the different particles is created in the folder `rundir/<exp_name>/output_fcost`. There is one fcost file per assimilation, with one line per particle.

If the program crashes or if you want more information, you can have a look at the log files, located in the folder `rundir/<exp_name>/output_log`. There is one log file per particle.

# References

- Dubinkina, S., Goosse, H., Sallaz-Damaz, Y., Crespin, E., and Crucifix, M.: Testing a particle filter to reconstruct climate changes over the past centuries. International Journal of Bifurcation and Chaos, 21, 3611–3618, [10.1142/S0218127411030763](10.1142/S0218127411030763), 2011.
- Dalaiden, Q., Goosse, H., Rezsöhazy, J., & Thomas, E. R. : Reconstructing atmospheric circulation and sea-ice extent in the west Antarctic over the past 200 years using data assimilation. Climate Dynamics, 57(11), 3479–3503. [https://doi.org/10.1007/s00382-021-05879-6](https://doi.org/10.1007/s00382-021-05879-6), 2021.

# Contributors

- François Klein
- Quentin Dalaiden