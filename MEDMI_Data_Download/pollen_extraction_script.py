from medmi_database import Dataset


Dates = ['2016-1-1 0','2019-12-31 23']
Dates_string = '2016-2019'


pollen_dict = {'Time range':Dates,'Latitude range': [48,60], 'Longitude range': [-11,3]}

dataset_base = 'midas.pollen_drnl_ob.{}'
file_base = 'data_pollen/pollen_{}_{}.csv'

pollen_species = ['alnus','ambrosia','artemisia','betula',\
					'corylus','fraxinus','platanus','poaceae',\
					'quercus','salix','ulmus','urtica']


print('extracting pollen data for dates: {} to {}'.format(Dates[0],Dates[1]))

for pspc in pollen_species:
	print('working on pollen species: {}'.format(pspc))
	pollen_dict['Source reference'] = dataset_base.format(pspc)
	fname = file_base.format(Dates_string,pspc)

	pollendata = Dataset(pollen_dict)
	pollendata.default()

	with open(fname,'w') as dfile:

		dfile.write('{} pollen daily count for date range: {} to {}\n'.format(pspc,Dates[0],Dates[1]))
		dfile.write('date,siteID,{}\n'.format(pspc))
	
		for data in pollendata.values():
			d_date = data['Time']
			d_siteid = data['Site identifier']	
			d_val = data['Value']
			dfile.write('{}, {}, {}\n'.format(d_date,d_siteid,d_val))


print('finished extracting pollen data')



