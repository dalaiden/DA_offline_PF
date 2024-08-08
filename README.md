# Offline Data Assimilation (particle filter; UCLouvain-ELIC)

Here is the code of the off-line data assimilation program based on the particle filter described, in
its online version, in Dubinkina et al. (2011). This program has been used in Klein et al. (2018, https://doi.org/10.1007/s00382-017-3853-0) and in Klein et al. (2019, https://doi.org/10.5194/cp-15-661-2019). Please see those publications for more information.  

It has been tested for several configurations but may need some adjustements depending on your system. Please contact François klein if you experience problems. 

The branch "master" only contains the scripts. The branch "daps" contains an example with all necessary inputs, that you should be able to run it directly once downloaded. This example corresponds to the experiment performed in the framework of the data assimation intercomparison project, initiated by the [DAPS](http://www.pages-igbp.org/1500-word-essay-pages.aspx) group with the support of PAGES.  

## Installation

- git clone https://gogs.elic.ucl.ac.be/TECLIM/DA_Offline.git

From the folder `DA_Offline`:
- Open the file `src/modules.load` and check that you have all requested libraries at the right place (currently templates are available for the clusters `storm`, `coriolis` and `cyclone`).
- Execute the Makefile: `cd src; . modules.load; make`.
- Check that you have your executables: `src/moddata_co` and `src/PartFilter`.

## Test experiment
Once you have compiled the sources, you can run a test experiment by executing the script `assim_offline.sc` (in the `DA_Offline` folder):

`./assim_offline.sc`

This example corresponds to the experiment performed in the framework of the data assimation intercomparison project, initiated by the [DAPS](http://www.pages-igbp.org/1500-word-essay-pages.aspx) group with the support of PAGES. Surface temperature pseudodata is assimilated over the period 1900-2000 using a LOVECLIM model ensemble.  

Information about the progress of the program will appear on the screen. Once finished, output files containing the fcosts of the different particles are created in the folder `rundir/<exp_name>/output_fcost`.

:warning: **The following lines need to be updated** :warning: They still contain useful information but correspond to a former version of DA_Offline.

## Experiment preparation

### Input files

Before giving where these input files have to be noted down and their dimension, here is a list of the input files absolutely necessary when launching a new experiment, and of the ones optional. There a three files required:
- The model ensemble member(s).
- The file containing the data that should be assimilated.
- The file containing the error of the data assimilated.

In addition, the following files can optionally be specified:
- If you want the particle filter to work with anomalies, two more files are needed, containing the references for computing the anomalies for the data and the ensemble members.
- A matrix of covariance for the variable assimilated can also be used.
- A file containing the domains if spatial averages are needed, if those domains are not rectangular.
- A file containing the standardized areas of each grid cells if the grid is irregular and if spatial averaged are performed.

These files included their full path must be specified in namelist_brut, except the model ensemble simulations (included in assim_offline.sc, see Section 2.2). All the optional variables can be left empty. For normal use, jump directly to line 27 of
namelist_brut:

- **OceanAreaPath** (optional) is the first variable potentially to fill in. OceanAreaPath should contain the netcdf file with the full path containing the variable areacello, which is the area of each grid cell of the file, standardized between 0 and 1. The
dimension of the file is lat*lon.

- **Box file** (optional) should point to a 2-dimensional (lat*lon) netcdf file (with full path) with, for each model grid cell, contains a number corresponding to the domain over which the average has to be made. If no spatial average is needed or if the domains are rectangular, leave this variable empty.

- **Line 44-61** Template showing the order in which the information about the variables assimilated should be defined between the lines 65 to 80. One bloc per variable assimilated should be filled in, with always the same structure:

    - var name (mandatory) is the name of the variable to be assimilated, as specified in the netcdf files.
    - var type (mandatory) is the realm of the variable assimilated, ie. atmos or ocean.
    - 4D structure (mandatory) is the dimension of the netcdf file. Do not change, it should always be time*lat*lon.
    - weight (mandatory) is the weight that is given, when more than one variable is assimilated, to each variable assimilated. Set weight to 1 if only one variable is assimilated.
    - sqrt(Ci) (mandatory) is the error of the data if it is uniform. Pay attention, this option has not been tested yet. Currently, the error are provided externally in a netcdf file (see below ErrObs).
    - Sigma (mandatory) Weight for the model covariance matrix. 0 if cov matrix should be diagonal. 
    - Number of domains (mandatory) Number of domains used. Set to 1 if no spatial averages are performed.
    - Box Type (optional) If your domain is rectangular, Box Type is set to square. If your box is irregular, Box Type is set to file, whose address is specified in line 28 in box_file.
    - ObsFileStart: startDateY (mandatory) Must be 1.
    - ObsFileStart: startDateD (mandatory) Must be 1.
    - DataObs (mandatory) Netcdf file including the data that should be assimilated, with full path. Dimension should be time*lat*lon.
    - RefObs (optional) Netcdf file of the reference file of the data assimilated, with full path. Usually, it is the average of every January values, February values, ... December values. The size of this file is 12*lat*lon.
    - ErrObs (mandatory) Netcdf file containing the error of the data assimilated, with full path. The dimension is lat*lon.
    - RefMod (optional) The reference file of the model simulations, with full path. It is usually the average of every January values, February values, ... December values for every ensemble members. The size of this file is 12*lat*lon.
    - Datacov (optional) Netcdf file containing the matrix of covariance for the variable assimilated.

- boxvar (optional) When using domains that are rectangular (when Box Type is set to square), this last variable includes the indices corresponding to each domain. The first numbers are the coordinates of the grid points where the average as to be written (lat, lon), and the last four correspond to the four corners of the box.

### Experimental design

The timing of the experiment and the options related to the ensemble size are specified in the file assim_offline.sc. For normal use, you are not supposed to change anything beyond the line 30. Here are the parameters to edit:

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

Once the first assimilation is finished, an output file containing the the fcosts of the different particles is created in the folder `rundir/<exp_name>/output_fcost`. There is one fcost file per assimilation, with one line per particle. A basic fcost
looks like:

|NAME | START | END |     FCOST     | CURRENT WEIGHT | TARGET WEIGHT | RESTART COPY | NEW WEIGHT|
| --- | ----- | --- | ------------- | -------------- | ------------- | ------------ | --------- |
| 1   |   1   |  12 | 0.9792574E-01 |       1        |      2        |      1       |      1    |
| 2   |  13   |  24 | 0.2758664E+00 |       1        |      0        |      2       |      1    |
| 3   |  25   |  36 | 0.7439515E-02 |       1        |      0        |      1       |      1    |
| 4   |  37   |  48 | 0.6619669E-02 |       1        |      0        |      2       |      1    |
| 5   |  49   |  60 | 0.2795728E+00 |       1        |      2        |      5       |      1    |
| 6   |  61   |  72 | 0.2261626E-01 |       1        |      0        |      5       |      1    |
| 7   |  73   |  84 | 0.7401583E+00 |       1        |      5        |      7       |      1    |
| 8   |  85   |  96 | 0.5692443E-01 |       1        |      0        |      8       |      1    |
| 9   |  97   | 108 | 0.3835564E-02 |       1        |      0        |      7       |      1    |

- **NAME** Number to distinguish the ensemble members, from 1 to the ensemble size.
- **START** Start of the period based on which the fcosts are computed.
- **END** End of the period based on which the fcosts are computed.
- **FCOST** Likelihood of the particle.
- **CURRENT WEIGHT** Always 1, not used here.
- **TARGET WEIGHT** The weight of the particle, needed to compute the weighted mean leading to the final reconstruction. In this example, the reconstruction for the assimilation step is (part1*2 + part5*2 + part7*5 ) / 9.
- **RESTART** COPY Only necessary in the online version.
- **NEW WEIGHT** Only necessary in the online version.

If the program crashes or if you want more information, you can have a look at the log files, located in the folder `rundir/<exp_name>/output_log`. There is one log file per particle.

## References
Dubinkina, S., Goosse, H., Sallaz-Damaz, Y., Crespin, E., and Crucifix, M.: Testing a particle filter to reconstruct climate changes over the past centuries, International Journal of Bifurcation and Chaos, 21, 3611–3618, doi:10.1142/S0218127411030763,2011.