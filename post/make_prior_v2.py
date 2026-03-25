#python3

import os, sys, errno
import numpy as np
from netCDF4 import Dataset
import shutil

#-------------
# Parameters |
#-------------
exp_id = str(sys.argv[-1])

loc_data_prior = '/cyfast/hxue/Data_processed/LEs' 
# loc_data_prior = '/cyfast/dalaiden/LEs/processed'

fname_model_ID = '../rundir/{}/info_prior/prior'.format(exp_id)
model_ID = open(fname_model_ID, 'r').read()

fname_year_a_prior = '../rundir/{}/info_prior/year_a_prior'.format(exp_id)
year_a_prior = int(open(fname_year_a_prior, 'r').read())

fname_year_b_prior = '../rundir/{}/info_prior/year_b_prior'.format(exp_id)
year_b_prior = int(open(fname_year_b_prior, 'r').read())
year_a_ano_prior = np.copy(year_a_prior)
year_b_ano_prior = np.copy(year_b_prior)
var_list = {
	'PRECT'                         : { 'var_ID'  : 'PRECT', 
										'unit_s'  : 'm/s'},
	'PSL'                           : { 'var_ID'  : 'PSL', 
										'unit_s'  : 'Pa'},
	'TREFHT'                        : { 'var_ID'  : 'TREFHT', 
										'unit_s'  : 'K'},
	'VPD'                         : { 'var_ID'  : 'VPD', 
										'unit_s'  : 'hpa'}
	'PDSI'                         : { 'var_ID'  : 'PDSI', 
										'unit_s'  : 'unitless'}										
}
list_seasons = ['JJA','MAMJJAS','ANN']
# season_id = 'MAMJJAS'
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
	nb_members=100
	years_prior = np.arange(1850, 2099+1)
elif model_ID == 'NorCPM1':
	nb_members=30
	years_prior = np.arange(1850, 2014+1)

print('create the prior for {} over {}-{}'.format(model_ID, year_a_prior, year_b_prior))

# Loop on all the variables to create
for var in var_list:

	print('	- '+var)

	var_ID = var_list[var]['var_ID']
	unit_s = var_list[var]['unit_s']

	# Loop on seasons
	for i_season in range(len(list_seasons)):
		
		season_id = list_seasons[i_season]

		print('   - '+season_id)

		# Load data
		fname = '{}/{}/{}/{}_{}-LE_{}_{}-{}.nc'.format(loc_data_prior, model_ID, var_ID, var_ID, model_ID, season_id, years_prior[0], years_prior[-1])
		nc = Dataset(fname)
		data_all = nc.variables[var][:]
		lat_raw, lon_raw = nc.variables['lat'][:], nc.variables['lon'][:]
		lon, lat = np.meshgrid(lon_raw, lat_raw)
		nc.close()

		data_all[np.abs(data_all) > 100000000] = np.nan


		for imember in range(data_all.shape[0]):
			
			data_m = data_all[imember,...].squeeze()

			# Compute the anomalies
			data_m = data_m - np.nanmean(data_m[(years_prior >= year_a_ano_prior) & (years_prior <= year_b_ano_prior),...], axis=0)

			# Keep the year_a - year_b period
			data_m = data_m[(years_prior >= year_a_prior) & (years_prior <= year_b_prior),...]

			# Export to netcdf
			output_dir = '{}/{}'.format(pid_s, var_ID)
			if (imember == 0) & (i_season == 0):
				if os.path.exists(output_dir):
					shutil.rmtree(output_dir)
			mkdir_p(output_dir)
			outfname = '{}/{}_prior_r{}i1p1_{}-{}_{}.nc'.format(output_dir, var_ID, str(imember+1).zfill(3), year_a_prior, year_b_prior, season_id)
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
