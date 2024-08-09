#python3

import os, sys, errno
import numpy as np
from netCDF4 import Dataset
import shutil

#-------------
# Parameters |
#-------------
exp_id = str(sys.argv[-1])

loc_data_prior = '/Users/dalaiden/Documents/DA_offline_PF/data_DA_offline_PF/prior'

fname_model_ID = '../rundir/{}/info_prior/prior'.format(exp_id)
model_ID = open(fname_model_ID, 'r').read()

fname_year_a_prior = '../rundir/{}/info_prior/year_a_prior'.format(exp_id)
year_a_prior = int(open(fname_year_a_prior, 'r').read())

fname_year_b_prior = '../rundir/{}/info_prior/year_b_prior'.format(exp_id)
year_b_prior = int(open(fname_year_b_prior, 'r').read())
year_a_ano_prior = np.copy(year_a_prior)
year_b_ano_prior = np.copy(year_b_prior)
var_list = {
	'PSL'                           : { 'var_ID'  : 'PSL', 
										'unit_s'  : 'Pa'},
	'TREFHT'                        : { 'var_ID'  : 'TREFHT', 
										'unit_s'  : 'K'},
	'ICEFRAC'                       : { 'var_ID'  : 'ICEFRAC', 
										'unit_s'  : 'ratio (0->)'},
	'sea-ice-extent_regions_RH2014' : { 'var_ID'  : 'sea-ice-extent_regions_RH2014', 
										'unit_s'  : '10^6 km^2'},
	'SAM_diff'                      : { 'var_ID'  : 'SAM_diff', 
										'unit_s'  : 'unitless'},
	'SST'                           : { 'var_ID'  : 'SST', 
										'unit_s'  : 'K'},
}
#-------------

pid_s = os.getpid()

f = open('pid_s', 'w')
f.write(str(pid_s))
f.close()

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc:  # Python ≥ 2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		# possibly handle other errno cases here, otherwise finally:
		else:
			raise

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
	nb_members=99
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


print('create the prior for {} over {}-{}'.format(model_ID, year_a_prior, year_b_prior))

# Loop on all the variables to create
for var in var_list:

	print('	- '+var)

	var_ID = var_list[var]['var_ID']
	unit_s = var_list[var]['unit_s']

	# Load data
	fname = '{}/{}/{}/{}_{}-LE_ANN_{}-{}.nc'.format(loc_data_prior, model_ID, var_ID, var_ID, model_ID, years_prior[0], years_prior[-1])
	nc = Dataset(fname)
	data_all = nc.variables[var][:]
	if (var_ID != 'sea-ice-extent_regions_RH2014') & (var_ID != 'sea-ice-area_regions_RH2014') & (var_ID != 'SAM_diff'):
		lat_raw, lon_raw = nc.variables['lat'][:], nc.variables['lon'][:]
		lon, lat = np.meshgrid(lon_raw, lat_raw)
	nc.close()

	data_all[np.abs(data_all) > 100000000] = np.nan

	if (var_ID == 'sea-ice-extent_regions_RH2014') | (var_ID == 'sea-ice-area_regions_RH2014'):

		data_all = np.concatenate((data_all, (data_all[:,:,0] + data_all[:,:,4])[:,:,None]), axis=-1) # West Antarctica
		data_all = np.concatenate((data_all, (data_all[:,:,1] + data_all[:,:,2] + data_all[:,:,3])[:,:,None]), axis=-1) # East Antarctica 
		data_all = np.concatenate((data_all, (data_all[:,:,2] + data_all[:,:,3])[:,:,None]), axis=-1) # 'King Hakon' + 'East Antarctica'

	for imember in range(data_all.shape[0]):
		
		data_m = data_all[imember,...].squeeze()

		# Compute the anomalies
		data_m = data_m - np.nanmean(data_m[(years_prior >= year_a_ano_prior) & (years_prior <= year_b_ano_prior),...], axis=0)

		# Keep the year_a - year_b period
		data_m = data_m[(years_prior >= year_a_prior) & (years_prior <= year_b_prior),...]

		# Export to netcdf
		output_dir = '{}/{}'.format(pid_s, var_ID)
		if imember == 0:
			if os.path.exists(output_dir):
				shutil.rmtree(output_dir)
		mkdir_p(output_dir)
		outfname = '{}/{}_prior_r{}i1p1_{}-{}_ANN.nc'.format(output_dir, var_ID, str(imember+1).zfill(3), year_a_prior, year_b_prior)
		# print('		export to: '+outfname)
		ncid = Dataset(outfname, 'w', format='NETCDF4')
		
		# define dimensions
		if len(np.shape(data_m)) == 2:
			lat_dim = ncid.createDimension('regions',data_m.shape[-1])
		elif len(np.shape(data_m)) == 3:
			dimid_lon = ncid.createDimension('lon', lon.shape[1])
			dimid_lat = ncid.createDimension('lat', lon.shape[0])
		dimid_time = ncid.createDimension('time', int(year_b_prior-year_a_prior+1))

		# define variables
		if len(np.shape(data_m)) == 2:
			lat_attr = ncid.createVariable('regions',np.float64,('regions',))
		elif len(np.shape(data_m)) == 3:
			varid_lon = ncid.createVariable('lon', 'f4', ('lon',))
			varid_lat = ncid.createVariable('lat', 'f4', ('lat',))
		varid_time = ncid.createVariable('time', 'f4', ('time',))
		if len(np.shape(data_m)) == 1:
			varid_zg = ncid.createVariable(var_ID, 'f4', ('time'))
		elif len(np.shape(data_m)) == 2:
			varid_zg = ncid.createVariable(var_ID, 'f4', ('time', 'regions'))
		elif len(np.shape(data_m)) == 3:
			varid_zg = ncid.createVariable(var_ID, 'f4', ('time', 'lat', 'lon'))

		if len(np.shape(data_m)) == 3:

			varid_lon.long_name = 'longitude coordinate'
			varid_lon.standard_name = 'longitude'
			varid_lon.units = 'degrees_east'
			varid_lon.axis = 'X'

			varid_lat.long_name = 'latitude coordinate'
			varid_lat.standard_name = 'latitude'
			varid_lat.units = 'degrees_north'
			varid_lat.axis = 'Y'

		varid_time.long_name = 'time'
		varid_time.standard_name = 'time'
		varid_time.units = 'months since 0001-01-01 00:00:00'
		varid_time.axis = '01-JAN-0000 00:00:00'

		varid_zg.long_name = var_ID
		varid_zg.standard_name = var_ID
		varid_zg.units = unit_s
		varid_zg.year_a_ano = year_a_ano_prior
		varid_zg.year_b_ano = year_b_ano_prior

		# write data
		if len(np.shape(data_m)) == 3:
			varid_lon[:] = lon[0, :]
			varid_lat[:] = lat[:, 0]

		varid_time[:] = np.arange(year_a_prior, year_b_prior+1)
		varid_zg[:] = data_m

		ncid.close()
