#Python3

import os, sys
import numpy as np
from netCDF4 import Dataset

"""

Several priors are available:

    1. CanESM2           1950 - 2100 (n = 50)
    2. CESM1             1920 - 2080 (n = 35)
    3. CESM1_LM          0850 - 2005 (n = 12)
    4. CESM2             1850 - 2100 (n = 100)
    5. CNRM-CM6-1        1850 - 2014 (n = 19)
    6. CSIRO-Mk3-6-0     1850 - 2100 (n = 30)
    7. GFDL-CM3          1920 - 2100 (n = 20)
    8. GFDL-ESM2M        1950 - 2100 (n = 30)
    9. IPSL-CM6A-LR      1850 - 2014 (n = 33)
   10. MPI-ESM           1850 - 2099 (n = 100)
   11. NorCPM1           1850 - 2014 (n = 30)
   12. CanESM5           1850 - 2014 (n = 40)
   13. MIROC6            1850 - 2014 (n = 50)
   14. UKESM1-0-LL       1850 - 2014 (n = 14)
   15. ACCESS-ESM1-5     1850 - 2014 (n = 40)

"""

# Parameters
var_list = {
    'var_tas_mid_lat_45NN' : { 'var_ID'  : 'temperature', 
                              'region_s': 'mid-latitudes',
                              'filter_s': 'Ryan_45NN'},
    'var_tas_ant_45NN'        : { 'var_ID'  : 'temperature', 
                              'region_s': 'ant',
                              'filter_s': 'filter_45NN'},
    'var_psl_mid_lat_45NN' : { 'var_ID'  : 'msl_pressure', 
                              'region_s': 'mid-latitudes',
                              'filter_s': 'Ryan_45NN'},
    'var_psl_ant_45NN'        : { 'var_ID'  : 'msl_pressure', 
                              'region_s': 'ant',
                              'filter_s': 'filter_45NN'},
    
}
# list_seasons = ["DJF", "MAM", "JJA", "SON"]
list_seasons = ["ANN"]
ocean_mask = False
loc_data = 'data'
loc_data_prior = '/Users/dalaiden/Documents/DA_offline_PF/data_DA_offline_PF/obs'
grid_size = '500km'
only_continuous_series = False
year_a_ano = 1958 # Obervations
year_b_ano = 2020 # Obervations
model_ID = 'CESM1_LM'
year_a_prior = 850
year_b_prior = 2005
year_a_ano_prior = 1958
year_b_ano_prior = 2005
interpolation_method = 'homemade'
error_inflation = np.arange(0.1, 10.1, 0.1)
default_value_constant_error = 1
#-----------------------

#-------------------------------------
# Export information about the prior |
#-------------------------------------
os.system('mkdir -p info_prior')

f = open('info_prior/prior', 'w')
f.write(model_ID)
f.close()

f = open('info_prior/year_a_prior', 'w')
f.write('{}'.format(str(year_a_prior)))
f.close()

f = open('info_prior/year_b_prior', 'w')
f.write('{}'.format(str(year_b_prior)))
f.close()

os.system('mkdir -p input')
os.chdir('input')

if model_ID == 'iCESM1':
    nb_members = 3
    years_prior = np.arange(851,2005+1)
elif model_ID == 'CanESM2':
    nb_members=50
    years_prior = np.arange(1950, 2100+1)
elif model_ID == 'CESM1':
    nb_members=35
    years_prior = np.arange(1850, 2100+1)
elif model_ID == 'CESM1_LM':
    nb_members=12
    years_prior = np.arange(850, 2005+1)
elif model_ID == 'CESM2':
    nb_members=100
    years_prior = np.arange(1850, 2100+1)
elif model_ID == 'CNRM-CM6-1':
    nb_members=19
    years_prior = np.arange(1850, 2014+1)
elif model_ID == 'CSIRO-Mk3-6-0':
    nb_members=30
    years_prior = np.arange(1850, 2100+1)
elif model_ID == 'GFDL-CM3':
    nb_members=20
    years_prior = np.arange(1920, 2100+1)
elif model_ID == 'GFDL-ESM2M':
    nb_members=30
    years_prior = np.arange(1950, 2100+1)
elif model_ID == 'IPSL-CM6A-LR':
    nb_members=33
    years_prior = np.arange(1850, 2014+1)
elif model_ID == 'MPI-ESM':
    nb_members=100
    years_prior = np.arange(1850, 2099+1)
elif model_ID == 'NorCPM1':
    nb_members=30
    years_prior = np.arange(1850, 2014+1)
elif model_ID == 'CanESM5':
    nb_members=40
    years_prior = np.arange(1850, 2014+1)
elif model_ID == 'MIROC6':
    nb_members=50
    years_prior = np.arange(1850, 2014+1)
elif model_ID == 'UKESM1-0-LL':
    nb_members=14
    years_prior = np.arange(1850, 2014+1)
elif model_ID == 'ACCESS-ESM1-5':
    nb_members=40
    years_prior = np.arange(1850, 2014+1)

# Loop on all the variables to create
for dir_out in var_list:

    var_ID = var_list[dir_out]['var_ID']
    region_s = var_list[dir_out]['region_s']
    filter_s = var_list[dir_out]['filter_s']

    print('Create files for {}'.format(dir_out))
    print(' - '+var_ID)
    print(' - '+region_s)
    print(' - '+filter_s)

    # Clean
    os.system('rm -rf {}'.format(dir_out))
    os.system('mkdir -p {}/data/files'.format(dir_out))
    os.system('mkdir -p {}/model/files/ensemble'.format(dir_out))

    #------
    # OBS |
    #------

    # Load the season(s)
    for i_season in range(1,len(list_seasons)+1):
        season_ID = list_seasons[i_season-1]

        fname = '{}/{}/{}/composites_{}_{}_filter_{}_1900-2025_grid_{}_{}_oceanic_grid_masked_{}_only_continuous_series_{}.nc'.format(loc_data, var_ID, grid_size, var_ID, region_s, filter_s, grid_size, season_ID, ocean_mask, only_continuous_series)

        nc = Dataset(fname)
        if i_season == 1:
            data_grid_tmp = nc.variables[var_ID][:]
            lon, lat = nc.variables['lon'][:], nc.variables['lat'][:]
        else:
            data_grid_tmp = nc.variables[var_ID][:]

        # Create the array with the four seasons
        if i_season == 1:
            data_grid = np.empty((len(list_seasons) * data_grid_tmp.shape[0], data_grid_tmp.shape[1], data_grid_tmp.shape[2])) * np.nan

        data_grid[i_season-1::len(list_seasons),:,:] = np.copy(data_grid_tmp)

        del data_grid_tmp

    # Reshape data
    data_rshp = np.empty((data_grid.shape[0], data_grid.shape[1] * data_grid.shape[2])) * np.nan
    lat_rshp = np.empty((data_grid.shape[1] * data_grid.shape[2])) * np.nan
    lon_rshp = np.empty((data_grid.shape[1] * data_grid.shape[2])) * np.nan

    k = 0
    for i_grid in range(data_grid.shape[1]):
        for j_grid in range(data_grid.shape[2]):
            data_rshp[:, k] = data_grid[:,i_grid,j_grid]
            lon_rshp[k] = lon[i_grid,j_grid]
            lat_rshp[k] = lat[i_grid,j_grid]
            k += 1

    # Create the mask to keep only the cell grid with data and not NaN
    mask = ~np.isnan(np.nanmean(data_rshp, axis=0))

    # # Save the mask
    # mask_obs = np.save('mask_obs.npy', mask)

    # Apply the mask on data
    matrice_results = data_rshp[:,mask]
    lon = lon_rshp[mask]
    lat = lat_rshp[mask]

    # Fix
    lon = np.arange(len(lon))
    lat = np.arange(len(lat))

    # Create new nc
    if len(list_seasons) == 4:
        outfile_name = '{}/data/files/{}_{}_1900-2025_4seasons.nc'.format(dir_out, var_ID, region_s)
    elif len(list_seasons) == 1:
        outfile_name = '{}/data/files/{}_{}_1900-2025.nc'.format(dir_out, var_ID, region_s)

    ncid = Dataset(outfile_name, 'w', format='NETCDF4')

    # Define dimensions
    dimid_dimsup = ncid.createDimension('dimsup', 1)
    dimid_lat = ncid.createDimension('lat', matrice_results.shape[1])
    dimid_lon = ncid.createDimension('lon', matrice_results.shape[1])
    dimid_time = ncid.createDimension('time', matrice_results.shape[0])

    varid_lat = ncid.createVariable('lat', 'f4', ('lat',))
    varid_lat.long_name = 'latitude coordinate'
    varid_lat.standard_name = 'latitude'
    varid_lat.units = 'degrees_north'
    varid_lat.axis = 'Y'
    varid_lat[:] = lat

    varid_lon = ncid.createVariable('lon', 'f4', ('lon',))
    varid_lon.long_name = 'longitude coordinate'
    varid_lon.standard_name = 'longitude'
    varid_lon.units = 'degrees_east'
    varid_lon.axis = 'X'
    varid_lon[:] = lon

    varid_time = ncid.createVariable('time', 'f4', ('time',))
    varid_time.long_name = 'time'
    varid_time.standard_name = 'time'
    varid_time.units = 'months since 0001-01-01 00:00:00'
    varid_time.axis = '01-JAN-0001 00:00:00'
    varid_time[:] = np.arange(1, matrice_results.shape[0] + 1)

    varid_d18O = ncid.createVariable('{}_{}'.format(var_ID, region_s), 'f4', ('time', 'lat', 'dimsup'))
    varid_d18O.long_name = '{}_{}'.format(var_ID, region_s)
    varid_d18O.standard_name = 'Obs {} {}'.format(var_ID, region_s)
    varid_d18O.units = ''
    varid_d18O.missing_value = -99.99
    matrice_results[np.isnan(matrice_results)] = -99.99
    varid_d18O[:] = matrice_results[:,:,None]

    ncid.close()

    print('     Netcdf with observations created')

    # Compute the reference for computing anomalies later
    nc = Dataset(outfile_name)
    data = nc.variables['{}_{}'.format(var_ID, region_s)][:]
    lon, lat = nc.variables['lon'][:], nc.variables['lat'][:]
    nc.close()
    data[(data >= -99.99001) & (data <= -99.98999)] = np.nan

    # Compute the mean over the specific period for each season
    year_tot = np.arange(1900, 2025+1)
    data_ref = np.empty((len(list_seasons), data.shape[1])) *np.nan
    for i_season in range(len(list_seasons)):
        data_season_tmp = data[i_season::len(list_seasons),:]
        data_ref[i_season, :] = np.nanmean(data_season_tmp[(year_tot >= year_a_ano) & (year_tot <= year_b_ano),:], axis=0).squeeze()

    matrice_results = data_ref[:,:,None]

    # create the nc
    if len(list_seasons) == 4:
        outfile_name = '{}/data/files/{}_{}_REF_4seasons.nc'.format(dir_out, var_ID, region_s)
    elif len(list_seasons) == 1:
        outfile_name = '{}/data/files/{}_{}_REF.nc'.format(dir_out, var_ID, region_s)
    ncid = Dataset(outfile_name, 'w', format='NETCDF4')

    # Define dimensions
    dimid_dimsup = ncid.createDimension('dimsup', 1)
    dimid_lat = ncid.createDimension('lat', matrice_results.shape[1])
    dimid_lon = ncid.createDimension('lon', matrice_results.shape[1])
    dimid_time = ncid.createDimension('time', matrice_results.shape[0])

    varid_lat = ncid.createVariable('lat', 'f4', ('lat',))
    varid_lat.long_name = 'latitude coordinate'
    varid_lat.standard_name = 'latitude'
    varid_lat.units = 'degrees_north'
    varid_lat.axis = 'Y'
    varid_lat[:] = lat

    varid_lon = ncid.createVariable('lon', 'f4', ('lon',))
    varid_lon.long_name = 'longitude coordinate'
    varid_lon.standard_name = 'longitude'
    varid_lon.units = 'degrees_east'
    varid_lon.axis = 'X'
    varid_lon[:] = lon

    varid_time = ncid.createVariable('time', 'f4', ('time',))
    varid_time.long_name = 'time'
    varid_time.standard_name = 'time'
    varid_time.units = 'months since 0001-01-01 00:00:00'
    varid_time.axis = '01-JAN-0001 00:00:00'
    varid_time[:] = np.arange(1, matrice_results.shape[0] + 1)

    varid_d18O = ncid.createVariable('{}_{}'.format(var_ID, region_s), 'f4', ('time', 'lat', 'dimsup'))
    varid_d18O.long_name = '{}_{}'.format(var_ID, region_s)
    varid_d18O.standard_name = 'Obs {}_{}'.format(var_ID, region_s)
    varid_d18O.units = ''
    varid_d18O.missing_value = -99.99
    varid_d18O[:] = matrice_results

    ncid.close()

    print('     Netcdf with the reference (for computing anomalies) created')

    # Create the netcdf containing the DA errors (constant error)
    for i_error in range(len(error_inflation)):

        if len(list_seasons) == 4:
            fname = '{}/data/files/{}_{}_REF_4seasons.nc'.format(dir_out, var_ID, region_s)
        elif len(list_seasons) == 1:
            fname = '{}/data/files/{}_{}_REF.nc'.format(dir_out, var_ID, region_s)
        nc = Dataset(fname)
        matrice_results = nc.variables['{}_{}'.format(var_ID, region_s)][:].squeeze()
        nc.close()

        if len(list_seasons) > 1:
            matrice_results = np.ones_like(matrice_results[0,...].squeeze()) * default_value_constant_error * error_inflation[i_error]
        else:
            matrice_results = np.ones_like(matrice_results.squeeze()) * default_value_constant_error * error_inflation[i_error]

        matrice_results = matrice_results[:,None]
        
        error_type = f'constant-error_factor-{error_inflation[i_error]:.2f}'.replace('.', '')
        if len(list_seasons) == 4:
            outfile_name = '{}/data/files/{}_{}_error_{}_4seasons.nc'.format(dir_out, var_ID, error_type, region_s)
        elif len(list_seasons) == 1:
            outfile_name = '{}/data/files/{}_{}_error_{}.nc'.format(dir_out, var_ID, error_type, region_s)
        ncid = Dataset(outfile_name, 'w', format='NETCDF4')

        dimid_dimsup = ncid.createDimension('dimsup', 1)
        dimid_lat = ncid.createDimension('lat', matrice_results.shape[0])
        dimid_lon = ncid.createDimension('lon', matrice_results.shape[0])

        varid_lat = ncid.createVariable('lat', 'f4', ('lat',))
        varid_lat.long_name = 'latitude coordinate'
        varid_lat.standard_name = 'latitude'
        varid_lat.units = 'degrees_north'
        varid_lat.axis = 'Y'
        varid_lat[:] = lat

        varid_lon = ncid.createVariable('lon', 'f4', ('lon',))
        varid_lon.long_name = 'longitude coordinate'
        varid_lon.standard_name = 'longitude'
        varid_lon.units = 'degrees_east'
        varid_lon.axis = 'X'
        varid_lon[:] = lon

        varid_rms = ncid.createVariable('rms', 'f4', ('lat', 'dimsup'))
        varid_rms.long_name = 'rms'
        varid_rms.standard_name = 'rms'
        varid_rms.units = ''
        varid_rms.missing_value = -99.99
        varid_rms[:] = matrice_results

        ncid.close()

        del matrice_results

    # Create the netcdf containing the DA errors (variance error)
    for i_error in range(len(error_inflation)):

        # Load observations
        if len(list_seasons) == 4:
            fname = '{}/data/files/{}_{}_1900-2025_4seasons.nc'.format(dir_out, var_ID, region_s)
        elif len(list_seasons) == 1:
            fname = '{}/data/files/{}_{}_1900-2025.nc'.format(dir_out, var_ID, region_s)

        nc = Dataset(fname)
        matrice_results = nc.variables['{}_{}'.format(var_ID, region_s)][:].squeeze()
        nc.close()

        matrice_results[(matrice_results >= -99.99001) & (matrice_results <= -99.98999)] = np.nan

        # Compute the std over 1981-2010
        if len(list_seasons) == 4:

            # Compute the std for each season
            for i_season in range(len(list_seasons)):
                data_tmp = matrice_results[i_season::len(list_seasons),...]
                std_season = np.nanstd(data_tmp[(year_tot >= 1958) & (year_tot <= 2023),:], axis=0)
                if i_season == 0:
                    std_all = np.ones_like(matrice_results)
                std_all[i_season::len(list_seasons),...] = std_all[i_season::len(list_seasons),...] * std_season * error_inflation[i_error]

            del matrice_results
            matrice_results = np.copy(std_all)
            del std_all

        elif len(list_seasons) == 1:
            matrice_results = np.nanstd(matrice_results[(year_tot >= 1958) & (year_tot <= 2023),:], axis=0) * error_inflation[i_error]
            matrice_results = matrice_results[:,None]
        
        error_type = f'std-error_factor-{error_inflation[i_error]:.2f}'.replace('.', '')
        if len(list_seasons) == 4:
            outfile_name = '{}/data/files/{}_{}_error_{}_4seasons.nc'.format(dir_out, var_ID, error_type, region_s)
        elif len(list_seasons) == 1:
            outfile_name = '{}/data/files/{}_{}_error_{}.nc'.format(dir_out, var_ID, error_type, region_s)
        ncid = Dataset(outfile_name, 'w', format='NETCDF4')

        dimid_dimsup = ncid.createDimension('dimsup', 1)
        if len(list_seasons) > 1:
            dimid_lat = ncid.createDimension('lat', matrice_results.shape[1])
            dimid_lon = ncid.createDimension('lon', matrice_results.shape[1])
            dimid_time = ncid.createDimension('time', matrice_results.shape[0])
        else:
            dimid_lat = ncid.createDimension('lat', matrice_results.shape[0])
            dimid_lon = ncid.createDimension('lon', matrice_results.shape[0])

        varid_lat = ncid.createVariable('lat', 'f4', ('lat',))
        varid_lat.long_name = 'latitude coordinate'
        varid_lat.standard_name = 'latitude'
        varid_lat.units = 'degrees_north'
        varid_lat.axis = 'Y'
        varid_lat[:] = lat

        varid_lon = ncid.createVariable('lon', 'f4', ('lon',))
        varid_lon.long_name = 'longitude coordinate'
        varid_lon.standard_name = 'longitude'
        varid_lon.units = 'degrees_east'
        varid_lon.axis = 'X'
        varid_lon[:] = lon

        if len(list_seasons) > 1:

            varid_time = ncid.createVariable('time', 'f4', ('time',))
            varid_time.long_name = 'time'
            varid_time.standard_name = 'time'
            varid_time.units = 'months since 0001-01-01 00:00:00'
            varid_time.axis = '01-JAN-0001 00:00:00'
            varid_time[:] = np.arange(1, matrice_results.shape[0] + 1)

        if len(list_seasons) > 1:

            varid_rms = ncid.createVariable('rms', 'f4', ('time','lat'))
            varid_rms.long_name = 'rms'
            varid_rms.standard_name = 'rms'
            varid_rms.units = ''
            varid_rms.missing_value = -99.99
            varid_rms[:] = matrice_results

        else:

            varid_rms = ncid.createVariable('rms', 'f4', ('lat', 'dimsup'))
            varid_rms.long_name = 'rms'
            varid_rms.standard_name = 'rms'
            varid_rms.units = ''
            varid_rms.missing_value = -99.99
            varid_rms[:] = matrice_results

        ncid.close()

        del matrice_results

    print('     All netcdf files containing errors exported')

    #--------
    # Prior |
    #--------

    if (interpolation_method == 'CDO') and (var_ID == 'temperature'):
        var_model = 'TREFHT'
    elif (interpolation_method == 'homemade') and (var_ID == 'temperature'):
        var_model = 'tas'
    elif (interpolation_method == 'homemade') and (var_ID == 'msl_pressure'):
        var_model = 'psl'

    if (ocean_mask == False) and (var_ID == 'temperature'):
        var_model = 'TREFHT'
    elif (ocean_mask == False) and (var_ID == 'msl_pressure'):
        var_model = 'PSL'

    # Loop on all the members
    for i_member in range(1, nb_members + 1):
        # Load the four seasons
        for i_season in range(len(list_seasons)):
            if (interpolation_method == 'homemade') and (ocean_mask == True):
                fname = '{}/{}/{}/{}_{}_{}_{}_{}_{}-{}_NEW.nc'.format(loc_data_prior, grid_size, model_ID, var_model, model_ID, str(i_member).zfill(3), grid_size, list_seasons[i_season], years_prior[0], years_prior[-1])
            else:
                if ocean_mask == False:
                    fname = '{}/{}/{}/{}_{}_{}_{}_{}_{}-{}_CDO_remapbil_with_ocean_data.nc'.format(loc_data_prior, grid_size, model_ID, var_model, model_ID, str(i_member).zfill(3), grid_size, list_seasons[i_season], years_prior[0], years_prior[-1])
                else:
                    fname = '{}/{}/{}/{}_{}_{}_{}_{}_{}-{}_{}.nc'.format(loc_data_prior, grid_size, model_ID, var_model, model_ID, str(i_member).zfill(3), grid_size, list_seasons[i_season], years_prior[0], years_prior[-1])

            # Load data
            nc = Dataset(fname, 'r')
            data_grid_tmp = nc.variables[var_model][:]
            nc.close()

            # Load geographical coordinates
            lon_lat_file = '/home/elic/dalaiden/20th_reconstruction_hgs_fogt/grids/grids/geographical_coordinates_SH_500km.nc'
            nc = Dataset(lon_lat_file, 'r')
            lon = nc.variables['lon'][:]
            lat = nc.variables['lat'][:]
            nc.close()

            data_grid_tmp[np.abs(data_grid_tmp) > 100000000] = np.nan

            if var_ID == 'msl_pressure':
                data_grid_tmp = data_grid_tmp / 100.0

            # Compute the anomalies here
            for i_season in range(len(list_seasons)):
                data_season_tmp = data_grid_tmp[i_season::len(list_seasons),...]
                ref = np.nanmean(data_season_tmp[(years_prior >= year_a_ano_prior) & (years_prior <= year_b_ano_prior),...], axis=0).squeeze()
                data_grid_tmp[i_season::len(list_seasons),...] = data_grid_tmp[i_season::len(list_seasons),...] - ref
                del data_season_tmp, ref

            # Keep the year_a_prior - year_b_prior period
            data_grid_tmp = data_grid_tmp[(years_prior >= year_a_prior) & (years_prior <= year_b_prior),...]

            # Create the array with the four seasons
            if i_season == 0:
                data_grid = np.empty((len(list_seasons) * data_grid_tmp.shape[0], data_grid_tmp.shape[1], data_grid_tmp.shape[2])) * np.nan

            data_grid[i_season::len(list_seasons),:,:] = np.copy(data_grid_tmp)
            del data_grid_tmp

        # Reshape data
        data_rshp = np.empty((data_grid.shape[0], data_grid.shape[1] * data_grid.shape[2])) * np.nan
        lat_rshp = np.empty((data_grid.shape[1] * data_grid.shape[2])) * np.nan
        lon_rshp = np.empty((data_grid.shape[1] * data_grid.shape[2])) * np.nan
        k = 0
        for i_grid in range(data_grid.shape[1]):
            for j_grid in range(data_grid.shape[2]):
                data_rshp[:, k] = data_grid[:,i_grid,j_grid]
                lon_rshp[k] = lon[i_grid,j_grid]
                lat_rshp[k] = lat[i_grid,j_grid]
                k += 1

        # Apply the mask
        matrice_results = data_rshp[:, mask]
        lat = lat_rshp[mask]
        lon = lon_rshp[mask]

        lon = np.arange(len(lon))
        lat = np.arange(len(lat))

        # Extract new nc's
        if len(list_seasons) == 4:
            outfile_name = '{}/model/files/ensemble/{}_{}_{}_{}_{}-{}_4seasons.nc'.format(dir_out, var_ID, region_s, model_ID, str(i_member).zfill(3), year_a_prior, year_b_prior)
        elif len(list_seasons) == 1:
            outfile_name = '{}/model/files/ensemble/{}_{}_{}_{}_{}-{}.nc'.format(dir_out, var_ID, region_s, model_ID, str(i_member).zfill(3), year_a_prior, year_b_prior)
        ncid = Dataset(outfile_name, 'w', format='NETCDF4')

        # Define dimensions
        dimid_dimsup = ncid.createDimension('dimsup', 1)
        dimid_lat = ncid.createDimension('lat', len(lat))
        dimid_lon = ncid.createDimension('lon', len(lon))
        dimid_time = ncid.createDimension('time', matrice_results.shape[0])

        varid_lat = ncid.createVariable('lat', 'f4', ('lat',))
        varid_lat.long_name = 'latitude coordinate'
        varid_lat.standard_name = 'ID point -> not real lat'
        varid_lat.units = 'degrees_north'
        varid_lat.axis = 'Y'
        varid_lat[:] = lat

        varid_lon = ncid.createVariable('lon', 'f4', ('lon',))
        varid_lon.long_name = 'longitude coordinate'
        varid_lon.standard_name = 'ID point -> not real lat'
        varid_lon.units = 'degrees_east'
        varid_lon.axis = 'X'
        varid_lon[:] = lon

        varid_time = ncid.createVariable('time', 'f4', ('time',))
        varid_time.long_name = 'time'
        varid_time.standard_name = 'time'
        varid_time.units = 'months since 0001-01-01 00:00:00'
        varid_time.axis = '01-JAN-0001 00:00:00'
        varid_time[:] = np.arange(1, matrice_results.shape[0] + 1)

        varid_d18O = ncid.createVariable('{}_{}'.format(var_ID, region_s), 'f4', ('time', 'lat', 'dimsup'))
        varid_d18O.long_name = '{}_{}'.format(var_ID, region_s)
        varid_d18O.standard_name = 'Simulated {}_{}'.format(var_ID, region_s)
        varid_d18O.units = ''
        varid_d18O.missing_value = -99.99
        varid_d18O[:] = matrice_results[:,:,None]

        ncid.close()

    print('     Netcdfs with simulated values created')

    # Compute the reference for computing anomalies later (NOT USED ANYMORE)
    year_tot = np.arange(year_a_prior,year_b_prior+1)
    for i_member in range(1, nb_members + 1):
        
        if len(list_seasons) == 4:
            fname = '{}/model/files/ensemble/{}_{}_{}_{}_{}-{}_4seasons.nc'.format(dir_out, var_ID, region_s, model_ID, str(i_member).zfill(3), year_a_prior, year_b_prior)
        elif len(list_seasons) == 1:
            fname = '{}/model/files/ensemble/{}_{}_{}_{}_{}-{}.nc'.format(dir_out, var_ID, region_s, model_ID, str(i_member).zfill(3), year_a_prior, year_b_prior)
        nc = Dataset(fname, 'r')
        data_d = nc.variables['{}_{}'.format(var_ID, region_s)][:]
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        nc.close()

        for i_season in range(len(list_seasons)):
            data_season_tmp = data_d[i_season::len(list_seasons),:]

            if (i_season == 0) and (i_member == 1):
                data_ref = np.empty((nb_members, len(list_seasons), data_season_tmp.shape[1])) * np.nan

            data_ref[i_member-1,i_season,:] = np.nanmean(data_season_tmp[(year_tot >= year_a_ano_prior) & (year_tot <= year_b_ano_prior), :], axis=0).squeeze()

    matrice_results = np.nanmean(data_ref, axis=0).squeeze() * 0

    # create the nc
    outfile_name = '{}/model/files/reference_model.nc'.format(dir_out)
    ncid = Dataset(outfile_name, 'w', format='NETCDF4')

    # Define dimensions
    dimid_dimsup = ncid.createDimension('dimsup', 1)
    if len(list_seasons) > 1:
        dimid_lat = ncid.createDimension('lat', matrice_results.shape[1])
        dimid_lon = ncid.createDimension('lon', matrice_results.shape[1])
        dimid_time = ncid.createDimension('time', matrice_results.shape[0])
    else:
        dimid_lat = ncid.createDimension('lat', matrice_results.shape[0])
        dimid_lon = ncid.createDimension('lon', matrice_results.shape[0])
        dimid_time = ncid.createDimension('time', 1)

    varid_lat = ncid.createVariable('lat', 'f4', ('lat',))
    varid_lat.long_name = 'latitude coordinate'
    varid_lat.standard_name = 'latitude'
    varid_lat.units = 'degrees_north'
    varid_lat.axis = 'Y'
    varid_lat[:] = lat

    varid_lon = ncid.createVariable('lon', 'f4', ('lon',))
    varid_lon.long_name = 'longitude coordinate'
    varid_lon.standard_name = 'longitude'
    varid_lon.units = 'degrees_east'
    varid_lon.axis = 'Y'
    varid_lon[:] = lon


    varid_time = ncid.createVariable('time', 'f4', ('time',))
    varid_time.long_name = 'time'
    varid_time.standard_name = 'time'
    varid_time.units = 'months since 0000-01-01 00:00:00'
    varid_time.axis = '01-JAN-0000 00:00:00'
    if len(list_seasons) > 1:
        varid_time[:] = np.arange(1, matrice_results.shape[0] + 1)
    else:
        varid_time[:] = 1

    
    varid_d18O = ncid.createVariable('{}_{}'.format(var_ID, region_s), 'f4', ('time', 'lat', 'dimsup'))
    varid_d18O.long_name = '{}_{}'.format(var_ID, region_s)
    varid_d18O.standard_name = 'Simulated {}_{}'.format(var_ID, region_s)
    varid_d18O.units = ''
    varid_d18O.missing_value = -99.99
    if len(list_seasons) > 1:
        varid_d18O[:] = matrice_results[:,:,None]
    else:
        varid_d18O[:] = matrice_results[None,:,None]

    ncid.close()

    print('     Netcdf with the reference (for computing anomalies) created')
