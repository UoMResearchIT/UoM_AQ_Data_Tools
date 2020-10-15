from medmi_database import Dataset
from cmath import polar
import argparse



def extraction_function(source_dict,settings_dict):
	datadata = Dataset(source_dict)
	datadata.default()
	
	with open(settings_dict['fname'],'w') as dfile:
	
		dfile.write(settings_dict['headstring'])
		dfile.write(settings_dict['columnstring'])
	
		for data in datadata.values():
			d_date = data['Time']
			d_siteid = data['Site identifier']	
			d_val = data['Value']
			dfile.write('{}, {}, {}\n'.format(d_date,d_siteid,d_val))


def extraction_add_data_function(source_dict,settings_dict,extra_datasets):
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
			d_extra = data['Additional field values']
			
			dfile.write('{}, {}, {}'.format(d_date,d_siteid,d_val))
			for d_ev in d_extra:
				dfile.write(', {}'.format(d_ev))
			dfile.write('\n')


def extraction_wind_function(source_dict,settings_dict):
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
	DEFAULT_DATE_RANGE = ['2016-1-1 0','2019-12-31 23']
	DEFAULT_DATES_STRING = '2016-2019'
	DEFAULT_LATITUDES = [48,60]
	DEFAULT_LONGITUDES = [-11,3]
	DEFAULT_VERBOSE = 0

	### general help text
	parser = argparse.ArgumentParser(
		description="*** A script for automated downloading of MET data for a given date range. ***")

	### parameters

	# filenames
	parser.add_argument("--outdir_prefix", "-o", dest="outdir_prefix",
						help="prefix to be appended on output directory name. Default: " + DEFAULT_OUT_DIR_PREFIX)
	parser.set_defaults(output_prefix=DEFAULT_OUT_DIR_PREFIX)

	# Dates
	parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',  help="start and end dates \
	 					(array - first two values only). Default: " + str(DEFAULT_DATE_RANGE))
	parser.add_argument("--dates_string", "-s", dest="dates_string", help= "date string to be used in file search. \
																		   Default: " + DEFAULT_DATES_STRING)
	parser.set_defaults(date_range=DEFAULT_DATE_RANGE)
	parser.set_defaults(dates_string=DEFAULT_DATES_STRING)

	# Latitude / longitude
	parser.add_argument("--latitude_range", "-t", dest="latitude_range", type=int, nargs='+',
						help="start and end latitude range (array - first two values only). Default: " + DEFAULT_LATITUDES)
	parser.set_defaults(latitude_range=DEFAULT_LATITUDES)
	parser.add_argument("--longitude_range", "-n", dest="longitude_range", type=int, nargs='+',
						help="start and end longitude range (array - first two values only). Default: " + DEFAULT_LONGITUDES)
	parser.set_defaults(longitude_range=DEFAULT_LONGITUDES)

	# Log verbose-ness
	parser.add_argument("--verbose", "-v", type=int,
						help="Level of output for debugging (Default: " + str(DEFAULT_VERBOSE) + " (0 = verbose output))")

	### Process Inputs
	args = parser.parse_args()

	if args.outdir_prefix:
		outdir_prefix = args.outdir_prefix
	else:
		print('No outdir_prefix given, so will use default: ' + DEFAULT_OUT_DIR_PREFIX)
		outdir_prefix = DEFAULT_OUT_DIR_PREFIX

	if args.date_range:
		date_range = args.date_range
	else:
		print('No date_range provided, so using default: ' + str(DEFAULT_DATE_RANGE))
		date_range = DEFAULT_DATE_RANGE

	if args.dates_string:
		dates_string = args.dates_string
	else:
		print('No dates_string provided, so using default: ' + str(DEFAULT_DATES_STRING))
		dates_string = DEFAULT_DATES_STRING

	if args.latitude_range:
		latitude_range = args.latitude_range
	else:
		print('No latitude_range provided, so using default: ' + str(DEFAULT_LATITUDES))
		latitude_range = DEFAULT_LATITUDES

	if args.longitude_range:
		longitude_range = args.longitude_range
	else:
		print('No longitude_range provided, so using default: ' + str(DEFAULT_LONGITUDES))
		longitude_range = DEFAULT_LONGITUDES

	if args.verbose:
		VERBOSE = max(args.verbose, 0)
		print('verbose: ', VERBOSE)
	else:
		print('No verbose flag provided, so using default: ' + str(DEFAULT_VERBOSE))
		VERBOSE = DEFAULT_VERBOSE



	dict_base = {'Time range':date_range,'Latitude range': [48,60], 'Longitude range': [-11,3]}
	
	
	print('extracting rain data for date range: {} to {}'.format(date_range[0],date_range[1]))
	rain_dict = dict_base.copy()
	rain_dict.update({'Source reference':'midas.rain_drnl_ob.prcp_amt 1'})
	rain_settings = {'fname':outdir_prefix + 'rain/rain_{}.csv'.format(dates_string),\
					'headstring':'Rain gauge daily data for date range: {} to {}\n'.format(date_range[0],date_range[1]),\
					'columnstring':'date,siteID,rain\n'}
	extraction_function(rain_dict,rain_settings)


	print('extracting temperature, relative humidity, pressure, and wet bulb temperature data for date range: {} to {}'
		  .format(date_range[0],date_range[1]))
	temperature_dict = dict_base.copy()
	temperature_dict.update({'Source reference':'midas.weather_hrly_ob.air_temperature'})
	temperature_settings = {'fname':outdir_prefix + 'met/temp_rh_press_wbulb_{}.csv'.format(dates_string),\
					'headstring':'Temperature, Relative Humidity, Station Pressure, and Wet Bulb Temperature hourly \
						data for date range: {} to {}\n'.format(date_range[0],date_range[1]),\
					'columnstring':'date,siteID,temperature,rh,pressure,wbtemp\n'}
	extra_datasets = ['midas.weather_hrly_ob.rltv_hum','midas.weather_hrly_ob.stn_pres','midas.weather_hrly_ob.dewpoint']
	extraction_add_data_function(temperature_dict,temperature_settings,extra_datasets)

	

	print('extracting wind data for date range: {} to {}'.format(date_range[0],date_range[1]))
	wind_dict = dict_base.copy()
	wind_dict.update({'Source reference':'midas.wind_mean_ob.mean_wind_speed 1','Complex wind type': True})
	wind_settings = {'fname':outdir_prefix + 'met/wind_{}.csv'.format(dates_string),\
					'headstring':'Wind speed and direction hourly data for date range: {} to {}\n'
						.format(date_range[0],date_range[1]),\
					'columnstring':'date,siteID,windspeed,winddir\n'}
	extraction_wind_function(wind_dict,wind_settings)


	print('all data extraction finished successfully')

