'''

Script for extracting data from the WRF meteorological files at given
measurement sites.

The WRF files are organised as 1 per day, with instantaneous output every 3 hours.

This requires the Salem library (installed via pip).

Variables to extract:
wind speed and direction (from U10, V10)



'''

import xarray as xr
import pandas as pd
from pathlib import Path
from numpy import datetime64 as dt

from wrf import getvar, ALL_TIMES, interpline, CoordPair, xy_to_ll, ll_to_xy 
from netCDF4 import Dataset


def load_sites_and_obtain_their_grid_locations(wrflist,sensors_file,site_list=[]):

    sens = pd.read_csv(sensors_file,usecols=['long','lat','sensor_name'])

    sens['xarray'],sens['yarray'] = ll_to_xy(wrflist,sens['lat'],sens['long'])
    # indexing later will be done:
    #    ds.T2C[time,yarray,xarray]

    # retain only those sites in our list (if the list isn't empty)
    if site_list:
        sens = sens[sens.sensor_name.isin(site_list)]


    # check to make sure that all stations are within our model domain - drop those that aren't

    if any(sens['xarray']>wrflist[0].dimensions['west_east'].size) or any(sens['yarray']>wrflist[0].dimensions['south_north'].size):
        print('dropping these stations outside of the model domain:')
        print(sens[sens['yarray']>wrflist[0].dimensions['south_north'].size])
        print(sens[sens['xarray']>wrflist[0].dimensions['west_east'].size])
        sens = sens[sens['yarray']<=wrflist[0].dimensions['south_north'].size]
        sens = sens[sens['xarray']<=wrflist[0].dimensions['west_east'].size]

    return(sens)


def create_dataframe_for_sites(wrflist,sens):


    # get the dates for the period of interest, using the keys from the resample function
    dates = getvar(wrflist,"times", timeidx=ALL_TIMES, method="cat").values


    # get the list of sites from the sensor data
    sensor_list = sens.sensor_name

    #%% create dataframe for storing model data
    #     this will start off with the format:
    #    rows: [date]
    # columns: [sensor name]

    model_data = pd.DataFrame(index=dates,columns=sensor_list)

    return(model_data)


def extract_wrf_data(wrflist,sens):


    # create the dataframe for our site data (species list is hardcoded in this function)
    windspeed_data = create_dataframe_for_sites(wrflist,sens)
    winddir_data = create_dataframe_for_sites(wrflist,sens)


    #%% run the function
    uvmet10_wspd_wdir = getvar(wrflist, "uvmet10_wspd_wdir", timeidx=ALL_TIMES, method="cat")


    for siteID in windspeed_data.columns:
        x_index = sens[sens.sensor_name==siteID].xarray.values
        y_index = sens[sens.sensor_name==siteID].yarray.values

        for t_index,dstamp in enumerate(windspeed_data.index):
            windspeed_data.loc[dstamp][siteID] = uvmet10_wspd_wdir[0,t_index,y_index,x_index].values.flatten()[0]
            winddir_data.loc[dstamp][siteID] = uvmet10_wspd_wdir[1,t_index,y_index,x_index].values.flatten()[0]


    return(windspeed_data,winddir_data)




if __name__ == '__main__':

    import argparse
    
    # read arguments from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--sensor_file","-s",help="path and name of the sensor locations file")
    parser.add_argument("--wrf_file_path","-wp",help="path to the WRF files")
    parser.add_argument("--date_range","-sd",help="year and month to process, YYYY-MM")
    parser.add_argument("--out_file_path","-op",help="path to the output directory")

    # read arguments from the command line
    args = parser.parse_args()


    if args.sensor_file:
        sensor_file = Path(args.sensor_file)
    else:
        #sensor_file = Path("./sensors_detailed.csv")
        print('Please provide sensor file location, using --sensor_file')

    if args.wrf_file_path:
        wrf_file_path = Path(args.wrf_file_path)
    else:
        print('Please provide wrf file path location, using --wrf_file_path')

    if args.date_range:
        date_range = args.date_range
    else:
        print('Please provide date range, using --date_range')

    if args.out_file_path:
        out_file_path = Path(args.out_file_path)
    else:
        print('Please provide output file path, using --out_file_path')


    site_list = ["Manchester Piccadilly [AQ]","Manchester Sharston [AQ]"]


    # operational code
    
    print('opening WRF files:')
    # create ordered list of files
    file_in_list = []
    if(wrf_file_path.is_dir()):
        for fil in wrf_file_path.glob('wrfout_d01_{}*'.format(date_range)): 
            file_in_list.append(fil) 
    else:
        print('{} is not a directory, no WRF files can be found')
    file_in_list.sort()
    
    print(file_in_list)
    # load files
    wrflist = []
    for fil in file_in_list:
        wrflist.append(Dataset(fil))
    
    print('loading sensor site locations')
    # load the sensor information 
    sens = load_sites_and_obtain_their_grid_locations(wrflist,sensor_file,site_list)

    print('extracting WRF data')
    # extract the data required from the WRF file
    (windspeed_dataframe,winddir_dataframe) = extract_wrf_data(wrflist,sens)

    print('writing to csv file')
    # write to file
    windspeed_dataframe.to_csv(out_file_path / 'wrf_windspeed_{}.csv'.format(date_range) )
    winddir_dataframe.to_csv(out_file_path / 'wrf_winddirection_{}.csv'.format(date_range) )
    

    print('finished')












