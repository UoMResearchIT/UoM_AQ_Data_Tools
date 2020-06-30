#
# Extract meteorological data for each of the pollen measurement sites
#

from medmi_database import Dataset
from cmath import polar


pollen_sites = {}
pollen_sites['BELFAST'] 		= (54.58537,   -5.93784, 15.0)
pollen_sites['BEVERLEY']        = (53.84243, -0.4252, 25.0)
pollen_sites['CARDIFF']         = (51.49536,   -3.21084, 25.0)
pollen_sites['CHESTER']         = (53.19835,   -2.89674, 28.0)
pollen_sites['ESKDALEMUIR']     = (55.31184, -3.20545, 236.0)
pollen_sites['EXETER']          = (50.73599, -3.53101, 10.0)
pollen_sites['IPSWICH']         = (52.05638, 1.1997, 41.0)
pollen_sites['LEICESTER']       = (52.62155,   -1.12227, 91.0)
pollen_sites['LONDON']          = (51.51022,   -0.11533, 45.0)
pollen_sites['MYLNEFIELD']      = (56.45699,   -3.07182, 31.0)
pollen_sites['PLYMOUTH']        = (50.37461,   -4.13725, 45.0)
pollen_sites['WIGHT']           = (50.71052,   -1.29944, 32.0)
pollen_sites['WORCESTER']       = (52.19716, -2.24165, 40.0)
pollen_sites['YORK']            = (53.94819,   -1.05194, 25.0)

Dates = ['2016-1-1 0','2019-12-31 23']
Dates_string = '2016-2019'


spatial_processing  = {'Method':'sp_idw_mean', 'Radius':50000}
spatial_temperature_processing  = {'Method':'sp_altadj_idw_mean', 'Radius':50000}
temporal_processing = {'Method':'tp_mean', 'Period':1} 

file_base = 'data_pollen/site_{}_{}_meteorology.txt'

# import the met datasets
rh = Dataset({'Source reference': 'midas.weather_hrly_ob.rltv_hum'})
rh.process(spatial_processing)
rh.process(temporal_processing)

temp = Dataset({'Source reference': 'midas.weather_hrly_ob.air_temperature'})
temp.process(spatial_temperature_processing)
temp.process(temporal_processing)

wind = Dataset({'Source reference': 'midas.wind_mean_ob.mean_wind_speed 1', \
				'Complex wind type': True})
wind.process(spatial_processing)
wind.process(temporal_processing)

rain = Dataset({'Source reference':'midas.rain_drnl_ob.prcp_amt 1'})
rain.process(spatial_processing)
rain.process(temporal_processing)


print('extracting averaged daily met data for pollen stations for date range: {} to {}'.format(Dates[0],Dates[1]))

for pollen_key in pollen_sites:
	print('working on site: {}'.format(pollen_key))
	site_lat = pollen_sites[pollen_key][0]
	site_lon = pollen_sites[pollen_key][1]
	ts = Dataset({'Source reference': 'MEDMI TS',\
				'Time range':Dates,\
				'Period':1,\
				'Latitude':site_lat,\
				'Longitude':site_lon,\
				'Altitude':pollen_sites[pollen_key][2]\
				})
				
	ts.link(rh)
	ts.link(temp)
	ts.link(wind)
	ts.link(rain)
	
	with open(file_base.format(pollen_key,Dates_string),'w') as dfile:
	
		dfile.write('{}, Latitude: {}, Longitude: {}\n'.format(pollen_key,site_lat,site_lon))
		dfile.write('date,relhum,temperature,windspeed,winddir,rain\n')
		
		for data in ts.values():
			d_date = data['Linked data'][0][0]['Time']
			d_rh = data['Linked data'][0][0]['Value']	
			d_temp = data['Linked data'][1][0]['Value']
			(d_wspd,d_wdir) = polar(data['Linked data'][2][0]['Value'])
			d_wdir = d_wdir*57.29577951308232
			if d_wdir < 0: d_wdir = d_wdir + 360.
			d_rain = data['Linked data'][3][0]['Value']
			dfile.write('{}, {}, {}, {}, {}, {}\n'.format(d_date,d_rh,d_temp,d_wspd,d_wdir,d_rain))
	
	
	
print('finished extracting met data for pollen stations')
