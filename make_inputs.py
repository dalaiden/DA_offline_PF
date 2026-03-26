#Python3

import os, sys
import numpy as np
from netCDF4 import Dataset
import pandas as pd
sys.path.insert(0, '/home/elic/dalaiden/python/')
sys.path.insert(0, '/elic/home/dalaiden/python/')
import func_Q as fq

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

"""

# Parameters
var_list = {
	'var_TRW_pr2'          : { 'var_ID'   : 'precipitation_Medi',
							  'obs_file' : 'supply_Medi_20240720.xlsx',
							  'folder_error' : 'results_20240720_prcp_1X1',
							  },
	'var_TRW_pr1'          : { 'var_ID'   : 'precipitation',
							  'obs_file' : 'Precipitation_records_20240619_1x1.xlsx',
							  'folder_error' : 'results_20240720_prcp_1X1',
							  },
	'var_TRW_tas'         : { 'var_ID'       : 'temperature', 
							  'obs_file'     : 'Temperature_records_20240619_1X1.xlsx',
							  'folder_error' : 'results_20240619_tem_1X1',
							  }					  
	
}

loc_data = '/home/elic/hxue/TRW_PSM/TRW_records'
loc_data_prior = '/cyfast/dalaiden/LEs/processed'
loc_PSM_results = '/home/elic/hxue/TRW_PSM'
year_a_ano = 1901 # Observations
year_b_ano = 2000 # Observations
model_ID = 'CESM1_LM'
year_a_prior = 850
year_b_prior = 1850
year_a_ano_prior = 850
year_b_ano_prior = 1850
error_inflation = np.arange(0.1, 10.1, 0.1)
default_value_constant_error = 1
#-----------------------

list_months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

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
os.chdir('input') # change the path to 'input/'

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

# Loop on all the variables to create
for dir_out in var_list:

	var_ID = var_list[dir_out]['var_ID']
	fname_obs = var_list[dir_out]['obs_file']
	folder_error = var_list[dir_out]['folder_error']

	print('Create files for {}'.format(dir_out))
	print(' - '+var_ID)
	print(' - file: '+fname_obs)

	# Clean
	os.system('rm -rf {}'.format(dir_out))
	os.system('mkdir -p {}/data/files'.format(dir_out))
	os.system('mkdir -p {}/model/files/ensemble'.format(dir_out))

	#------
	# OBS |
	#------

	# Load observations
	fname = '{}/{}'.format(loc_data, fname_obs)
	dfs = pd.read_excel(fname, sheet_name='Metadata')
	lat_records = dfs.grid_lat.values
	lon_records = dfs.grid_lon.values
	site_names = dfs.Number.values
	del dfs
	dfs = pd.read_excel(fname, sheet_name='Values', header=None)
	# years_data = dfs.values[:,0]
	# data_records = dfs.values[:,1:]	
	years_data = dfs.values[999:,0]
	data_records = dfs.values[999:,1:]
	del dfs
	obs_start_yr = int(years_data[0])
	obs_end_yr = int(years_data[-1])

	#----------------------------------------
	# Put the data on the grid of the prior |
	#----------------------------------------

	# Load the grid
	grid_info = '{}/{}/TREFHT/TREFHT_{}-LE_ANN_{}-{}.nc'.format(loc_data_prior, model_ID, model_ID, years_prior[0], years_prior[-1])
	nc = Dataset(grid_info)
	lon, lat = np.meshgrid(nc.variables['lon'][:], nc.variables['lat'][:])
	nc.close()

	# Load the land mask
	nc = Dataset('/cofast/dalaiden/20th_reconstruction_hgs_fogt/LEs/processed/fx/land_mask.nc')
	ocean_mask = nc.variables['mask'][:]
	nc.close()
	cont_mask = np.where(ocean_mask == 1, np.nan, 1)

	# Apply the mask
	lon, lat = lon * cont_mask, lat * cont_mask

	# Put on the data on the grid
	data_grid = np.empty((len(years_data), len(lon_records), lon.shape[0], lon.shape[1])) * np.nan

	for i_site in np.arange(len(lat_records)):

		# Compute distance
		dist_mat = np.empty_like(lat) * np.nan
		for i_lat in np.arange(lat.shape[0]):
			for i_lon in np.arange(lat.shape[1]):
				if np.isnan(lat[i_lat,i_lon]) == False:
					dist_mat[i_lat,i_lon] = fq.haversine((lat_records[i_site], lon_records[i_site]),(lat[i_lat,i_lon],lon[i_lat,i_lon]))

		# Find the min
		idx_lat, idx_lon = np.where(dist_mat == np.nanmin(dist_mat))
		idx_lat, idx_lon = idx_lat[0], idx_lon[0]

		# print('lat_record: {}; lat_grid: {}'.format(lat_records[i_site], lat[idx_lat, idx_lon]))
		# print('lon_record: {}; lon_grid: {}'.format(lon_records[i_site], lon[idx_lat, idx_lon]))

		# Store the data at the right location
		data_grid[:, i_site, idx_lat, idx_lon] = data_records[:,i_site].copy()

	# Average over the second dimension
	data_grid = np.nanmean(data_grid, axis=1)

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
	outfile_name = '{}/data/files/TRW_{}_{}-{}.nc'.format(dir_out, var_ID, obs_start_yr, obs_end_yr)

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

	varid_d18O = ncid.createVariable('TRW_{}'.format(var_ID), 'f4', ('time', 'lat', 'dimsup'))
	varid_d18O.long_name = 'TRW_{}'.format(var_ID)
	varid_d18O.standard_name = 'Obs TRW {}'.format(var_ID)
	varid_d18O.units = ''
	varid_d18O.missing_value = -99.99
	matrice_results[np.isnan(matrice_results)] = -99.99
	varid_d18O[:] = matrice_results[:,:,None]

	ncid.close()

	print('     Netcdf with observations created')

	# Compute the reference for computing anomalies later
	nc = Dataset(outfile_name)
	data = nc.variables['TRW_{}'.format(var_ID)][:]
	lon, lat = nc.variables['lon'][:], nc.variables['lat'][:]
	nc.close()
	data[(data >= -99.99001) & (data <= -99.98999)] = np.nan

	# Compute the mean over the specific period for each season
	year_tot = np.arange(obs_start_yr, obs_end_yr+1)
	data_ref = np.nanmean(data[(year_tot >= year_a_ano) & (year_tot <= year_b_ano),:], axis=0).squeeze()

	matrice_results = data_ref[None,:,None]

	# create the nc
	outfile_name = '{}/data/files/TRW_{}_REF.nc'.format(dir_out, var_ID)
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

	varid_d18O = ncid.createVariable('TRW_{}'.format(var_ID), 'f4', ('time', 'lat', 'dimsup'))
	varid_d18O.long_name = 'TRW_{}'.format(var_ID)
	varid_d18O.standard_name = 'Obs TRW_{}'.format(var_ID)
	varid_d18O.units = ''
	varid_d18O.missing_value = -99.99
	varid_d18O[:] = matrice_results

	ncid.close()

	print('     Netcdf with the reference (for computing anomalies) created')

	#---------------------------------|
	# Load the error based on the PSM |
	#---------------------------------|

	# Load the grid
	grid_info = '{}/{}/TREFHT/TREFHT_{}-LE_ANN_{}-{}.nc'.format(loc_data_prior, model_ID, model_ID, years_prior[0], years_prior[-1])
	nc = Dataset(grid_info)
	lon, lat = np.meshgrid(nc.variables['lon'][:], nc.variables['lat'][:])
	nc.close()

	# Load the land mask
	nc = Dataset('/cofast/dalaiden/20th_reconstruction_hgs_fogt/LEs/processed/fx/land_mask.nc')
	ocean_mask = nc.variables['mask'][:]
	nc.close()
	cont_mask = np.where(ocean_mask == 1, np.nan, 1)

	# Apply the mask
	lon, lat = lon * cont_mask, lat * cont_mask

	# Loop on sites on put the error on the grid of the model
	error_grid = np.empty_like(lon) * np.nan
	for i_site in range(len(site_names)):

		fname = '{}/{}/PSM_parameters_{}.txt'.format(loc_PSM_results, folder_error, site_names[i_site])
		TRW_model = pd.read_csv(fname, skiprows=2, header=None, sep=';')
		TRW_model = TRW_model.values

		error = TRW_model[-2,1]

		# Compute distance
		dist_mat = np.empty_like(lat) * np.nan
		for i_lat in np.arange(lat.shape[0]):
			for i_lon in np.arange(lat.shape[1]):
				if np.isnan(lat[i_lat,i_lon]) == False:
					dist_mat[i_lat,i_lon] = fq.haversine((lat_records[i_site], lon_records[i_site]),(lat[i_lat,i_lon],lon[i_lat,i_lon]))

		# Find the min
		idx_lat, idx_lon = np.where(dist_mat == np.nanmin(dist_mat))
		idx_lat, idx_lon = idx_lat[0], idx_lon[0]

		# print('lat_record: {}; lat_grid: {}'.format(lat_records[i_site], lat[idx_lat, idx_lon]))
		# print('lon_record: {}; lon_grid: {}'.format(lon_records[i_site], lon[idx_lat, idx_lon]))

		# Store the data at the right location
		error_grid[idx_lat, idx_lon] = np.copy(error)

	# Reshape data
	data_rshp = np.empty((error_grid.shape[0] * error_grid.shape[1])) * np.nan
	lat_rshp = np.empty((error_grid.shape[0] * error_grid.shape[1])) * np.nan
	lon_rshp = np.empty((error_grid.shape[0] * error_grid.shape[1])) * np.nan

	k = 0
	for i_grid in range(error_grid.shape[0]):
		for j_grid in range(error_grid.shape[1]):
			data_rshp[k] = error_grid[i_grid,j_grid]
			lon_rshp[k] = lon[i_grid,j_grid]
			lat_rshp[k] = lat[i_grid,j_grid]
			k += 1

	# Apply the mask on data
	matrice_results = data_rshp[mask]
	lon = lon_rshp[mask]
	lat = lat_rshp[mask]

	# Fix
	lon = np.arange(len(lon))
	lat = np.arange(len(lat))

	matrice_results = matrice_results[:,None]

	# Create the netcdf containing the DA errors (PSM error)
	for i_error in range(len(error_inflation)):
		
		error_type = f'PSM-error_factor-{error_inflation[i_error]:.2f}'.replace('.', '')
		outfile_name = '{}/data/files/TRW_{}_PSM_error_{}.nc'.format(dir_out, var_ID, error_type)
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
		varid_rms[:] = np.copy(matrice_results * error_inflation[i_error])

		ncid.close()

	print('     All netcdf files containing errors exported')


	#--------
	# Prior |
	#--------

	print(' - load the prior')

	# Workflow: load temperature and precipitation at monthly timescale but member by member

	# Loop on all the members
	for i_member in range(1, nb_members + 1):

		print('   - simu: {}'.format(i_member))

		print('     - temperature')
		
		#-------------------------------------|
		# Load temperature for all the months |
		#-------------------------------------|

		for i_month in range(len(list_months)):
			
			fname = '{}/{}/{}/TREFHT_{}-LE_{}_{}-{}.nc'.format(loc_data_prior, model_ID, 'TREFHT', model_ID, list_months[i_month], years_prior[0], years_prior[-1])

			# Load data
			nc = Dataset(fname)
			data_grid_tmp = nc.variables['TREFHT'][i_member-1,...].squeeze()
			if i_month == 0:
				lon, lat = np.meshgrid(nc.variables['lon'][:], nc.variables['lat'][:])
			nc.close()

			if i_month == 0:
				data_grid_all = np.empty((len(years_prior)*12, data_grid_tmp.shape[1], data_grid_tmp.shape[2])) * np.nan

			data_grid_all[i_month::12,...] = data_grid_tmp.copy()

			del data_grid_tmp

		# Now, we only want to keep temperature at record locations, so reshape
		data_rshp = np.empty((data_grid_all.shape[0], data_grid_all.shape[1] * data_grid_all.shape[2])) * np.nan
		lat_rshp = np.empty((data_grid_all.shape[1] * data_grid_all.shape[2])) * np.nan
		lon_rshp = np.empty((data_grid_all.shape[1] * data_grid_all.shape[2])) * np.nan

		k = 0
		for i_grid in range(data_grid_all.shape[1]):
			for j_grid in range(data_grid_all.shape[2]):
				data_rshp[:, k] = data_grid_all[:,i_grid,j_grid]
				lon_rshp[k] = lon[i_grid,j_grid]
				lat_rshp[k] = lat[i_grid,j_grid]
				k += 1

		# Apply the mask on data
		prior_tas_m_reshaped = data_rshp[:,mask]
		lon = lon_rshp[mask]
		lat = lat_rshp[mask]

		del data_grid_all # Clean

		#---------------------------------------|
		# Load precipitation for all the months |
		#---------------------------------------|

		print('     - precipitation')

		for i_month in range(len(list_months)):
			
			fname = '{}/{}/{}/PRECT_{}-LE_{}_{}-{}.nc'.format(loc_data_prior, model_ID, 'PRECT', model_ID, list_months[i_month], years_prior[0], years_prior[-1])

			# Load data
			nc = Dataset(fname)
			data_grid_tmp = nc.variables['PRECT'][i_member-1,...].squeeze()
			nc.close()

			if i_month == 0:
				data_grid_all = np.empty((len(years_prior)*12, data_grid_tmp.shape[1], data_grid_tmp.shape[2])) * np.nan

			data_grid_all[i_month::12,...] = data_grid_tmp.copy()

			del data_grid_tmp

		# Now, we only want to keep temperature at record locations, so reshape
		data_rshp = np.empty((data_grid_all.shape[0], data_grid_all.shape[1] * data_grid_all.shape[2])) * np.nan

		k = 0
		for i_grid in range(data_grid_all.shape[1]):
			for j_grid in range(data_grid_all.shape[2]):
				data_rshp[:, k] = data_grid_all[:,i_grid,j_grid]
				k += 1

		# Apply the mask on data
		prior_pr_m_reshaped = data_rshp[:,mask]

		del data_grid_all # Clean

		# ----------------------------------------------------
		print('       - load the results of the PSM to make the prior')
		# ----------------------------------------------------

		# Loop on all the sites
		for i_site in range(len(site_names)):

			# Find the name of the site
			lon_s, lat_s = lon[i_site],  lat[i_site]
			idx_site = np.where((lon_s == lon_records) & (lat_s == lat_records))

			site_name_s = site_names[idx_site[0]][0]

			# Load the results from the PSM
			fname = '{}/{}/PSM_parameters_{}.txt'.format(loc_PSM_results, folder_error, site_name_s)
			TRW_model = pd.read_csv(fname, skiprows=2, header=None, sep=';')
			TRW_model = TRW_model.values

			# Compute the modelled TRW based on PSM results
			nb_vars = TRW_model.shape[0]-2

			# Get the names of the variables
			TRW_prior = np.zeros_like(years_prior)
			for i_var_model in range(nb_vars):

				var_id = TRW_model[i_var_model,0].split('_')[-1]
				season_id = TRW_model[i_var_model,0].split('_')[0]

				if var_id == 'tas':

					# Compute the seasonal mean
					data_tas_reg = fq.annual_season(prior_tas_m_reshaped[:, i_site], season_id)

					# TRW component for this variable
					TRW_prior = TRW_prior + (data_tas_reg * TRW_model[i_var_model,-1])

				elif var_id == 'prcp':

					# Compute the seasonal mean
					data_prcp_reg = fq.annual_season(prior_pr_m_reshaped[:, i_site], season_id) * 1000

					# TRW component for this variable
					TRW_prior = TRW_prior + (data_prcp_reg * TRW_model[i_var_model,-1])

			# Store the data for each site
			if i_site == 0:
				matrice_results = np.empty((len(years_prior), len(lat_records))) * np.nan
			matrice_results[:,i_site] = np.copy(TRW_prior)

		# Compute anomalies
		matrice_results = matrice_results - np.mean(matrice_results[(years_prior >= year_a_ano_prior) & (years_prior <= year_b_ano_prior), :], axis=0)[None,:]

		# Keep the specific period
		matrice_results = matrice_results[(years_prior >= year_a_prior) & (years_prior <= year_b_prior),:]

		# Fix
		lon = np.arange(len(lon))
		lat = np.arange(len(lat))

		# Extract new nc's
		outfile_name = '{}/model/files/ensemble/TRW_{}_{}_{}_{}-{}.nc'.format(dir_out, var_ID, model_ID, str(i_member).zfill(3), year_a_prior, year_b_prior)
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

		varid_d18O = ncid.createVariable('TRW_{}'.format(var_ID), 'f4', ('time', 'lat', 'dimsup'))
		varid_d18O.long_name = 'TRW_{}'.format(var_ID)
		varid_d18O.standard_name = 'Simulated TRW_{}'.format(var_ID)
		varid_d18O.units = ''
		varid_d18O.missing_value = -99.99
		varid_d18O[:] = matrice_results[:,:,None]

		ncid.close()

	print('     Netcdfs with simulated values created')
	
	
	# Compute the reference for computing anomalies later (NOT USED ANYMORE)
	year_tot = np.arange(year_a_prior,year_b_prior+1)
	fname = '{}/model/files/ensemble/TRW_{}_{}_{}_{}-{}.nc'.format(dir_out, var_ID, model_ID, str(1).zfill(3), year_a_prior, year_b_prior)
	nc = Dataset(fname, 'r')
	data_d = nc.variables['TRW_{}'.format(var_ID)][:].squeeze()
	lat = nc.variables['lat'][:]
	lon = nc.variables['lon'][:]
	nc.close()

	matrice_results = np.nanmean(data_d, axis=0).squeeze() * 0

	# create the nc
	outfile_name = '{}/model/files/reference_model.nc'.format(dir_out)
	ncid = Dataset(outfile_name, 'w', format='NETCDF4')

	# Define dimensions
	dimid_dimsup = ncid.createDimension('dimsup', 1)
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
	varid_time[:] = 1


	varid_d18O = ncid.createVariable('TRW_{}'.format(var_ID), 'f4', ('time', 'lat', 'dimsup'))
	varid_d18O.long_name = 'TRW_{}'.format(var_ID)
	varid_d18O.standard_name = 'Simulated TRW_{}'.format(var_ID)
	varid_d18O.units = ''
	varid_d18O.missing_value = -99.99
	varid_d18O[:] = matrice_results[None,:,None]

	ncid.close()

	print('     Netcdf with the reference (for computing anomalies) created')
