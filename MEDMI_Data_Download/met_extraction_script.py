from medmi_database import Dataset
from cmath import polar
import argparse

import os, errno


def create_directory(dir_name):
	try:
		os.makedirs(dir_name)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise


def extract_data(source_dict, settings_dict, extra_datasets=[]):
	ADDITIONAL_VALS_KEY = 'Additional field values'
	datadata = Dataset(source_dict)
	datadata.default()
	for ds in extra_datasets:
		datadata.add(ds)
	
	with open(settings_dict['fname'],'w') as dfile:
		dfile.write(settings_dict['headstring'])
		dfile.write(settings_dict['columnstring'])
	
		for data in datadata.values():
			d_date = data['Time']
			d_siteid = data['Site identifier']	
			d_val = data['Value']

			if ADDITIONAL_VALS_KEY in data:
				d_extra = data[ADDITIONAL_VALS_KEY]
			else:
				d_extra = []
			
			dfile.write('{}, {}, {}'.format(d_date,d_siteid,d_val))
			for d_ev in d_extra:
				dfile.write(', {}'.format(d_ev))
			dfile.write('\n')


def extract_data_wind(source_dict,settings_dict):
	datadata = Dataset(source_dict)
	datadata.default()
	
	with open(settings_dict['fname'],'w') as dfile:
	
		dfile.write(settings_dict['headstring'])
		dfile.write(settings_dict['columnstring'])
	

		for data in datadata.values():
			d_date = data['Time']
			d_siteid = data['Site identifier']	
			(d_wspd,d_wdir) = polar(data['Value'])
			d_wdir = d_wdir*57.29577951308232
			if d_wdir < 0: d_wdir = d_wdir + 360.

			dfile.write('{}, {}, {}, {}\n'.format(d_date,d_siteid,d_wspd,d_wdir))




if __name__ == '__main__':
	global VERBOSE
	DEFAULT_OUT_DIR_PREFIX = 'data_'
	DEFAULT_OUT_FILE_SUFFIX = '2016-2019'
	DEFAULT_DATE_RANGE = ['2016-1-1 0','2019-12-31 23']
	DEFAULT_LATITUDES = [48,60]
	DEFAULT_LONGITUDES = [-11,3]
	DEFAULT_VERBOSE = 0

	### general help text
	parser = argparse.ArgumentParser(
		description="*** A script for automated downloading of MET data for a given date range. ***")

	### parameters

	# output directory/file names
	parser.add_argument("--outdir_prefix", "-o", dest="outdir_prefix", type=str,
						help="prefix to be appended on output directory name. Default: {}".format(DEFAULT_OUT_DIR_PREFIX))
	parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
						help="suffix to be appended to output file name. Default: {}".format(DEFAULT_OUT_FILE_SUFFIX))

	# Dates
	parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',  help="start and end dates \
	 					(array - first two values only). Default: {}".format(str(DEFAULT_DATE_RANGE)[1:-1].replace(',','')))


	# Latitude / longitude
	parser.add_argument("--latitude_range", "-t", dest="latitude_range", type=int, nargs='+',
						help="start and end latitude range (array - first two values only). \
							Default: {}".format(str(DEFAULT_LATITUDES)[1:-1].replace(',','')))
	parser.add_argument("--longitude_range", "-n", dest="longitude_range", type=int, nargs='+',
						help="start and end longitude range (array - first two values only). \
							Default: {}".format(str(DEFAULT_LONGITUDES)[1:-1].replace(',','')))

	# Log verbose-ness
	parser.add_argument("--verbose", "-v", type=int,
						help="Level of output for debugging (Default: {} (0 = verbose output))".format(str(DEFAULT_VERBOSE)))

	### Process Inputs
	args = parser.parse_args()

	if args.outdir_prefix:
		outdir_prefix = args.outdir_prefix
		print('Using outdir_prefix: {}'.format(outdir_prefix))
	else:
		print('No outdir_prefix given, so will use default: {}'.format(DEFAULT_OUT_DIR_PREFIX))
		outdir_prefix = DEFAULT_OUT_DIR_PREFIX

	if args.outfile_suffix:
		outfile_suffix = args.outfile_suffix
		print('Using outfile_suffix: {}'.format(outfile_suffix))
	else:
		print('No outfile_suffix provided, so using default: {}'.format(str(DEFAULT_OUT_FILE_SUFFIX)))
		outfile_suffix = DEFAULT_OUT_FILE_SUFFIX

	if args.date_range:
		if len(args.date_range) >= 2:
			date_range = args.date_range[0:2]
		else:
			raise ValueError(
				'Unable to obtain 2 dates from input --date_range: {}'.format(str(args.date_range)))
		print('Using date range: {}'.format(str(date_range)[1:-1].replace(',','')))
	else:
		print('No date_range provided, so using default: {}'.format(str(DEFAULT_DATE_RANGE)[1:-1].replace(',','')))
		date_range = DEFAULT_DATE_RANGE

	if args.latitude_range:
		if len(args.latitude_range) >= 2:
			latitude_range = args.latitude_range[0:2]
		else:
			raise ValueError('Unable to obtain 2 values from input --latitude_range: {}'.format(str(args.latitude_range)))
		print('Using latitude range: {}'.format(str(latitude_range)[1:-1].replace(',','')))
	else:
		print('No latitude_range provided, so using default: {}'.format(str(DEFAULT_LATITUDES)[1:-1].replace(',','')))
		latitude_range = DEFAULT_LATITUDES

	if args.longitude_range:
		if len(args.longitude_range) >= 2:
			longitude_range = args.longitude_range[0:2]
		else:
			raise ValueError('Unable to obtain 2 values from input --longitude_range: {}'.format(str(args.longitude_range)))
		print('Using longitude range: {}'.format(str(longitude_range)[1:-1].replace(',','')))
	else:
		print('No longitude_range provided, so using default: {}'.format(str(DEFAULT_LONGITUDES)[1:-1].replace(',','')))
		longitude_range = DEFAULT_LONGITUDES

	if args.verbose:
		VERBOSE = max(args.verbose, 0)
		print('verbose: ', VERBOSE)
	else:
		print('No verbose flag provided, so using default: {}'.format(str(DEFAULT_VERBOSE)))
		VERBOSE = DEFAULT_VERBOSE


	# Check directory name is valid and create directory
	new_dir = '{}met'.format(outdir_prefix)
	try:
		print('Creating directory: {}, unless it already exists.'.format(new_dir))
		create_directory(new_dir)
	except:
		raise ValueError('Unable to create directory/file: {}'.format(new_dir))


	### Prepare inputs and perform data extraction

	dict_base = {'Time range':date_range,'Latitude range': latitude_range, 'Longitude range': longitude_range}

	# Set file names
	filename_rain = '{}/rain_{}.csv'.format(new_dir, outfile_suffix)
	filename_temp_etc = '{}/temp_rh_press_wbulb_{}.csv'.format(new_dir, outfile_suffix)
	filename_wind = '{}/wind_{}.csv'.format(new_dir, outfile_suffix)
	filename_pollen = '{}/pollen_{}_{}.csv'.format(new_dir, outfile_suffix, '{}')

	# Rain
	print('extracting rain data for date range: {} to {}'.format(date_range[0],date_range[1]))
	rain_dict = dict_base.copy()
	rain_dict.update({'Source reference':'midas.rain_drnl_ob.prcp_amt 1'})
	rain_settings = {'fname':filename_rain,\
					'headstring':'Rain gauge daily data for date range: {} to {}\n'.format(date_range[0],date_range[1]),\
					'columnstring':'date,siteID,rain\n'}
	extract_data(rain_dict, rain_settings)

	# Temperature, humiidty, pressure, wet bulp temp
	print('extracting temperature, relative humidity, pressure, and wet bulb temperature data for date range: {} to {}'
		  .format(date_range[0],date_range[1]))
	temperature_dict = dict_base.copy()
	temperature_dict.update({'Source reference':'midas.weather_hrly_ob.air_temperature'})
	temperature_settings = {'fname':filename_temp_etc,\
					'headstring':'Temperature, Relative Humidity, Station Pressure, and Wet Bulb Temperature hourly \
						data for date range: {} to {}\n'.format(date_range[0],date_range[1]),\
					'columnstring':'date,siteID,temperature,rh,pressure,wbtemp\n'}
	extra_datasets = ['midas.weather_hrly_ob.rltv_hum','midas.weather_hrly_ob.stn_pres','midas.weather_hrly_ob.dewpoint']
	extract_data(temperature_dict, temperature_settings, extra_datasets)

	# Wind
	print('extracting wind data for date range: {} to {}'.format(date_range[0],date_range[1]))
	wind_dict = dict_base.copy()
	wind_dict.update({'Source reference':'midas.wind_mean_ob.mean_wind_speed 1','Complex wind type': True})
	wind_settings = {'fname':filename_wind,\
					'headstring':'Wind speed and direction hourly data for date range: {} to {}\n'
						.format(date_range[0],date_range[1]),\
					'columnstring':'date,siteID,windspeed,winddir\n'}
	extract_data_wind(wind_dict, wind_settings)

	# Pollen
	dataset_base = 'midas.pollen_drnl_ob.{}'
	pollen_species = ['alnus', 'ambrosia', 'artemisia', 'betula', 'corylus', 'fraxinus', 'platanus', 'poaceae',
					  'quercus', 'salix', 'ulmus', 'urtica']

	for pspc in pollen_species:
		print('working on pollen species: {}'.format(pspc))
		fname = filename_pollen.format(pspc)
		pollen_dict = dict_base.copy()
		pollen_dict['Source reference'] = dataset_base.format(pspc)
		pollen_settings = {'fname': fname,
						 'headstring': '{} pollen daily count for date range: {} to {}\n'.format(pspc, date_range[0], date_range[1]),
						 'columnstring': 'date,siteID,{}\n'.format(pspc)}
		extract_data(pollen_dict, pollen_settings)



