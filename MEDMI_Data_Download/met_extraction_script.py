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
	parser = argparse.ArgumentParser(
		description="*** A script for automated downloading of MET data for a given date range. ***")

	Dates = ['2016-1-1 0','2019-12-31 23']
	Dates_string = '2016-2019'
	dict_base = {'Time range':Dates,'Latitude range': [48,60], 'Longitude range': [-11,3]}
	
	
	print('extracting rain data for date range: {} to {}'.format(Dates[0],Dates[1]))
	rain_dict = dict_base.copy()
	rain_dict.update({'Source reference':'midas.rain_drnl_ob.prcp_amt 1'})
	rain_settings = {'fname':'data_rain/rain_{}.csv'.format(Dates_string),\
					'headstring':'Rain gauge daily data for date range: {} to {}\n'.format(Dates[0],Dates[1]),\
					'columnstring':'date,siteID,rain\n'}
	extraction_function(rain_dict,rain_settings)


	print('extracting temperature, relative humidity, pressure, and wet bulb temperature data for date range: {} to {}'.format(Dates[0],Dates[1]))
	temperature_dict = dict_base.copy()
	temperature_dict.update({'Source reference':'midas.weather_hrly_ob.air_temperature'})
	temperature_settings = {'fname':'data_met/temp_rh_press_wbulb_{}.csv'.format(Dates_string),\
					'headstring':'Temperature, Relative Humidity, Station Pressure, and Wet Bulb Temperature hourly data for date range: {} to {}\n'.format(Dates[0],Dates[1]),\
					'columnstring':'date,siteID,temperature,rh,pressure,wbtemp\n'}
	extra_datasets = ['midas.weather_hrly_ob.rltv_hum','midas.weather_hrly_ob.stn_pres','midas.weather_hrly_ob.dewpoint']
	extraction_add_data_function(temperature_dict,temperature_settings,extra_datasets)

	

	print('extracting wind data for date range: {} to {}'.format(Dates[0],Dates[1]))
	wind_dict = dict_base.copy()
	wind_dict.update({'Source reference':'midas.wind_mean_ob.mean_wind_speed 1','Complex wind type': True})
	wind_settings = {'fname':'data_met/wind_{}.csv'.format(Dates_string),\
					'headstring':'Wind speed and direction hourly data for date range: {} to {}\n'.format(Dates[0],Dates[1]),\
					'columnstring':'date,siteID,windspeed,winddir\n'}
	extraction_wind_function(wind_dict,wind_settings)


	print('all data extraction finished successfully')

