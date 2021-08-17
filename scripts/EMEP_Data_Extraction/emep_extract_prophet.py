'''

Script for extracting data from the EMEP data files at given
measurement sites.

The EMEP files are organised as 1 per 2 months, with instantaneous output every 1 hour.

This requires the Salem library (installed via pip).

Variables to extract:
NO2   - SURF_ppb_NO2
NOx   - SURF_ppb_NO2 + SURF_ppb_NO
SO2   - SURF_ppb_SO2
O3    - SURF_ppb_O3
PM2.5 - SURF_ug_PM25X_rh50
PM10  - SURF_ug_PM10


'''

import xarray as xr
import pandas as pd
from pathlib import Path
import salem
from numpy import datetime64 as dt


def load_sites_and_obtain_their_grid_locations(wrf_in,sensors_file,site_list):

    # open the sensor dataset
    sens = pd.read_csv(sensors_file,usecols=['long','lat','sensor_name'])

    # retain only those sites in our list
    sens = sens[sens.sensor_name.isin(site_list)]

    # get the indexes from the wrf file
    with salem.open_wrf_dataset(wrf_in) as ds:
        sens['iarray'],sens['jarray'] = ds.salem.grid.transform(sens['long'],sens['lat'],nearest=True)

    #%% check to make sure that all stations are within our model domain - drop those that aren't
    if any(sens['iarray']>ds.dims['west_east']) or any(sens['jarray']>ds.dims['south_north']):
        print('dropping these stations outside of the model domain:')
        print(sens[sens['jarray']>ds.dims['south_north']])
        print(sens[sens['iarray']>ds.dims['west_east']])
        sens = sens[sens['jarray']<=ds.dims['south_north']]
        sens = sens[sens['iarray']<=ds.dims['west_east']]

    return(sens)


def create_dataframe_for_sites(emep_ds,sens):

    # get the dates we are going to process (cutting off any dates for which we have less than 24 data points)
    dates = emep_ds.time.values     

    # get the list of sites from the sensor data
    sensor_list = sens.sensor_name

    ### create dataframe for storing model data
    #     this will start off with the format:
    #    rows: [date] 
    # columns: [sensor name]


    model_data = pd.DataFrame(index=dates,columns=sensor_list)

    return(model_data)



def extract_emep_data(emep_ds,sens):

    # create the dataframe for our site data (species list is hardcoded in this function)
    model_data = create_dataframe_for_sites(emep_ds,sens)

    # load the required data (hardcoded for the moment)
    NO2_data = emep_ds.SURF_ppb_NO2
    NOX_data = emep_ds.SURF_ppb_NO2 + emep_ds.SURF_ppb_NO



    # loop through sites, and dates of interest, to extract the required data
    for siteID in model_data.columns:
        i_index = sens[sens.sensor_name==siteID].iarray.values[0]
        j_index = sens[sens.sensor_name==siteID].jarray.values[0]


        for t_index,dstamp in enumerate(model_data.index):
            model_data.loc[dstamp][siteID] = NO2_data[t_index,j_index,i_index].values


    # drop the last row of the dataset, as this will be repeated in other datasets
    #   (this is a dodgy hardcoded solution...)
    model_data.drop(model_data.tail(1).index,inplace=True)  

    return(model_data)



if __name__ == '__main__':

    import argparse
    
    # read arguments from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--sensor_file","-s",help="path and name of the sensor locations file")
    parser.add_argument("--wrf_file","-w",help="path and name of the WRF output file for location information")
    parser.add_argument("--emep_file","-e",help="path and name of the EMEP output file to process")
    parser.add_argument("--out_file","-o",help="path and name of the EMEP output file to process")

    # read arguments from the command line
    args = parser.parse_args()

    if args.sensor_file:
        sensor_file = Path(args.sensor_file)
    else:
        #sensor_file = Path("./sensors_detailed.csv")
        print('Please provide sensor file location, using --sensor_file')

    if args.wrf_file:
        wrf_in = Path(args.wrf_file)
    else:
        #wrf_path = Path('../WRF_UK/2017/')
        #wrf_in = wrf_path / "wrfout_d01_2017-05-25_00:00:00"
        print('please provide WRF output file, using --wrf_file')

    if args.emep_file:
        emep_in = Path(args.emep_file)
    else:
        #emep_path = Path('../EMEP_UK/')
        #emep_in = emep_path / "Base_hourInst_MayJune2017.nc" 
        print('please provide EMEP output file, using --emep_file')

    if args.out_file:
        out_file = args.out_file
    else:
        #outpath = Path('./')
        #outfile = outpath / "emep_sitedata_MayJune2017.csv"
        print('using default output filename: ./outfile.csv')
        out_file = Path('./outfile.csv')

    site_list = ["Manchester Piccadilly [AQ]","Manchester Sharston [AQ]"]

    ### operational code

    print('loading sensor site locations')
    # open the sensor dataset, with grid locations
    sens = load_sites_and_obtain_their_grid_locations(wrf_in,sensor_file,site_list)

    print('opening EMEP file')
    # open EMEP file
    emep_ds = xr.open_dataset(emep_in)
    # indexing later will be done:
    #    emep_ds.T2C[time,jarray,iarray]  
    
    print('extracting EMEP data')
    # extract the model data for each site
    site_dataframe = extract_emep_data(emep_ds,sens)

    print('writing to csv file')
    # write the model data to csv
    site_dataframe.to_csv(out_file)

    print('finished')











