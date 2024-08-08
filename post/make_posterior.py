#python3

import numpy as np
from netCDF4 import Dataset
import sys, os, errno
sys.path.insert(0, '/home/elic/dalaiden/python/')
import func_Q as fq
from glob import glob
import time as time_pc
from joblib import Parallel, delayed
import shutil

# Parameters
nb_cores = 2 # 5
exp_id = str(sys.argv[-2])
outfolder = str(sys.argv[-1])
year_b = 2025
list_variables = [
				'PSL',
				'TREFHT',
				'SAM_diff',
				'sea-ice-extent_regions_RH2014',
				'sea-ice-area_regions_RH2014',
				'ICEFRAC',
				'PRECT',
				'SST',
				'U10m',
				'V10m',
				  ]
list_var_units = [
			  'hPa',
			  'K',
			  'unitless',
			  '10^6 km^2',
			  '10^6 km^2',
			  'ratio: 0->1',
			  'm/s',
			  'K',
			  'm s**-1',
			  'm s**-1',
				   ]

fname_model_ID = '../info_prior/prior'
model_ID = open(fname_model_ID, 'r').read()
#-----------------

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc:  # Python ≥ 2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		# possibly handle other errno cases here, otherwise finally:
		else:
			raise

# Load the pid
pid_s = int(open('pid_s', 'r').read())

# Load the first year of the reconstruction
year_a = int(open('../rundir/{}/first_year_rec'.format(exp_id), 'r').read())

print('Make the posterior ({})'.format(outfolder))

# Loop on all the variables
for i_var in range(len(list_variables)):

	var_ID = list_variables[i_var]
	var_unit = list_var_units[i_var]

	print('	- '+var_ID)

	# Create the folder
	mkdir_p(outfolder)
	mkdir_p('{}/{}'.format(outfolder, exp_id))

	season_ID = 'ANN'

	# Load all the simulations
	list_files = glob('{}/{}/*.nc'.format(pid_s, var_ID))
	list_files.sort()
	for i_file in range(len(list_files)):

		nc = Dataset(list_files[i_file])
		data_tmp = nc.variables[var_ID][:]
		if (i_file == 0) & (len(np.shape(data_tmp)) > 2):
			lat, lon = nc.variables['lat'][:], nc.variables['lon'][:]
			if len(lat.shape) > 2:
				print('There is a problem in the lat dimension')
		nc.close()

		if len(np.shape(data_tmp)) == 1:
			data_tmp = data_tmp[:,None,None]
		elif len(np.shape(data_tmp)) == 2:
			data_tmp = data_tmp[:,:,None]

		if i_file == 0:
			data = np.empty((len(list_files), data_tmp.shape[0], data_tmp.shape[1], data_tmp.shape[2])) * np.nan

		data[i_file,...] = np.copy(data_tmp)

		del data_tmp

	# Don't keep latitudes > 20deg South
	if 'lat' in locals():
		
		lat_init, lon_init = np.copy(lat), np.copy(lon)

		msk_lat = lat < -20

		data, lat = data[:,:,msk_lat,:], lat[msk_lat]

	# List all the fcosts
	list_fcost = glob('../rundir/{}/output_fcost/fcost_*'.format(exp_id))
	list_fcost = fq.natural_sort(list_fcost)
	nb_years = len(list_fcost)

	# Function to reconstruct the year
	def reconstruct_year(i_time):

		# Load the file
		fid = open(list_fcost[i_time], 'r')
		fcost_output = [line.strip().split() for line in fid.readlines()[1:]]
		fid.close()

		part_weight = np.array([float(x[5]) for x in fcost_output])
		timing_start = np.array([int(x[1]) for x in fcost_output]) - 1
		timing_end = np.array([int(x[2]) for x in fcost_output])

		particles_kept_rel = np.sum(part_weight!=0) / len(part_weight)

		series = np.tile(np.arange(len(list_files)), int(part_weight.shape[0]/len(list_files)))

		# Compute the weighted mean + weighted std
		nb_particles = len(part_weight)
		weighted_sum = np.zeros((data.shape[2], data.shape[3]))
		idx_std = 0
		weighted_std_tmp = np.empty((nb_particles, data.shape[2], data.shape[3])) * np.nan
		for i_part in range(nb_particles):

			weighted_sum = weighted_sum + part_weight[i_part] * np.nanmean(data[series[i_part],timing_start[i_part]:timing_end[i_part],...], axis=0)

			# Weighted std
			if part_weight[i_part] > 0:
				for i_std in range(int(part_weight[i_part])):
					
					weighted_std_tmp[idx_std,...] = np.nanmean(data[series[i_part],timing_start[i_part]:timing_end[i_part],...], axis=0)

					idx_std += 1

		# Now compute the final weigthted mean and std
		weighted_mean_time = (weighted_sum/nb_particles).squeeze()
		weighted_std_time = np.nanstd(weighted_std_tmp,axis=0).squeeze()

		return weighted_mean_time, weighted_std_time, particles_kept_rel

	# Run it
	results = Parallel(n_jobs=nb_cores)(delayed(reconstruct_year)(i_time) for i_time in  fq.progressbar(range(nb_years)))

	weighted_mean = np.empty((nb_years, data.shape[2], data.shape[3])) * np.nan
	weighted_std = np.empty((nb_years, data.shape[2], data.shape[3])) * np.nan
	particles_kept_rel_all = np.empty((nb_years)) * np.nan
	for i_time in range(nb_years):
		if len(np.shape(results[i_time][0])) == 1:
			weighted_mean[i_time,...] = results[i_time][0][:,None]
			weighted_std[i_time,...] = results[i_time][1][:,None]
		elif len(np.shape(results[i_time][0])) == 0:
			weighted_mean[i_time,...] = results[i_time][0][None,None]
			weighted_std[i_time,...] = results[i_time][1][None,None]
		elif len(np.shape(results[i_time][0])) == 2:
			weighted_mean[i_time,...] = results[i_time][0]
			weighted_std[i_time,...] = results[i_time][1]
		else:
			print('issue...')
			sys.exit()
		particles_kept_rel_all[i_time] = results[i_time][2]

	del results

	# Squeeze
	weighted_std = weighted_std.squeeze()
	weighted_mean = weighted_mean.squeeze()

	# Get back to the full domain
	if 'lat' in locals():

		weighted_mean_new = np.empty((weighted_mean.shape[0],lat_init.size,lon_init.size)) * np.nan
		weighted_std_new = np.empty((weighted_std.shape[0],lat_init.size,lon_init.size)) * np.nan

		weighted_mean_new[:,msk_lat,:] = np.copy(weighted_mean)
		weighted_std_new[:,msk_lat,:] = np.copy(weighted_std)

		lat = np.copy(lat_init)
		weighted_mean = np.copy(weighted_mean_new)
		weighted_std = np.copy(weighted_std_new)

		del lat_init, weighted_mean_new, weighted_std_new

	# Export the file
	outfile_name = '{}/{}/{}_DA-{}_{}-{}_{}.nc'.format(outfolder, exp_id, var_ID, model_ID, year_a, year_b, season_ID)
	# print('		--> '+outfile_name)
	outfile = Dataset(outfile_name, mode='w',format='NETCDF4')
	outfile.history = 'Created on '+time_pc.ctime(time_pc.time())
	outfile.program_used = 'Python'
	outfile.source = os.path.abspath(__file__)
	outfile.author = 'Q. Dalaiden (quentin.dalaiden@uclouvain.be)'

	# # Create the dimensions
	time_dim  = outfile.createDimension('time',None)
	if len(np.shape(weighted_std)) == 3:
		lat_dim = outfile.createDimension('lat',lat.size)
		lon_dim = outfile.createDimension('lon',lon.size)
	elif len(np.shape(weighted_std)) == 2:
		reg_dim = outfile.createDimension('regions',weighted_mean.shape[1])

	time_attr               = outfile.createVariable('time',np.float64,('time',))
	time_attr.standard_name = 'time'
	time_attr.long_name     = 'time'
	time_attr.units         = 'years since {}-01-01'.format(year_a)
	time_attr.calendar      = 'standard'
	time_attr.axis          = 'T'

	if len(np.shape(weighted_std)) == 3:

		lat_attr           = outfile.createVariable('lat',np.float64,('lat'))
		lat_attr.long_name = 'latitude'
		lat_attr.units     = 'degrees_north'

		lon_attr           = outfile.createVariable('lon',np.float64,('lon'))
		lon_attr.long_name = 'longitude'
		lon_attr.units     = 'degrees_east'

	if len(np.shape(weighted_std)) == 3:

		temp_var_nc  = outfile.createVariable('weighted_mean',
			                                      datatype='f4',
			                                      dimensions=('time','lat','lon'),
			                                      fill_value=1e+30)                                      
		temp_var_nc.long_name = var_ID
		temp_var_nc.units     = var_unit

		temp_var_nc2  = outfile.createVariable('weighted_std',
			                                      datatype='f4',
			                                      dimensions=('time','lat','lon'),
			                                      fill_value=1e+30)                                      
		temp_var_nc2.long_name = var_ID
		temp_var_nc2.units     = var_unit

	elif len(np.shape(weighted_std)) == 2:

		temp_var_nc  = outfile.createVariable('weighted_mean',
			                                      datatype='f4',
			                                      dimensions=('time','regions'),
			                                      fill_value=1e+30)                                      
		temp_var_nc.long_name = var_ID
		temp_var_nc.units     = var_unit

		temp_var_nc2  = outfile.createVariable('weighted_std',
			                                      datatype='f4',
			                                      dimensions=('time','regions'),
			                                      fill_value=1e+30)                                      
		temp_var_nc2.long_name = var_ID
		temp_var_nc2.units     = var_unit

	elif len(np.shape(weighted_std)) == 1:

		temp_var_nc  = outfile.createVariable('weighted_mean',
			                                      datatype='f4',
			                                      dimensions=('time'),
			                                      fill_value=1e+30)                                      
		temp_var_nc.long_name = var_ID
		temp_var_nc.units     = var_unit

		temp_var_nc2  = outfile.createVariable('weighted_std',
			                                      datatype='f4',
			                                      dimensions=('time'),
			                                      fill_value=1e+30)                                      
		temp_var_nc2.long_name = var_ID
		temp_var_nc2.units     = var_unit


	temp_var_nc3  = outfile.createVariable('particles_kept_rel',
		                                      datatype='f4',
		                                      dimensions=('time'))                                      
	temp_var_nc3.long_name = 'Relative number of particles kept '
	temp_var_nc3.units     = '0->1'

	# Fill the nc
	time_attr[:] = np.arange(nb_years)
	if len(np.shape(weighted_std)) == 3:
		lon_attr[:]  = lon
		lat_attr[:]  = lat
	weighted_mean[np.isnan(weighted_mean)] = 1e+30
	weighted_std[np.isnan(weighted_std)] = 1e+30
	temp_var_nc[:] = weighted_mean
	temp_var_nc2[:] = weighted_std
	temp_var_nc3[:] = particles_kept_rel_all

	outfile.close()

	if 'lat' in locals():
		del lat, lon

# Remove the prior folder
shutil.rmtree(str(pid_s))
