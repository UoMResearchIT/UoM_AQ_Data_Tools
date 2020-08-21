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
import pyreadr
from pathlib import Path
import salem
from numpy import datetime64 as dt


def load_sites_and_obtain_their_grid_locations(wrf_in,sensors_file,sensor_file_type):

    # open the sensor dataset
    if(sensor_file_type == 'CSV'):
        sens = pd.read_csv(sensors_file,usecols=['long','lat','sensor_name'])
    elif(sensor_file_type == 'RDATA'):
        metadata = pyreadr.read_r(sensors_file.as_posix())
        sens = metadata['AURN_metadata'][['site_id','latitude','longitude']].drop_duplicates() 
        sens = sens.rename(columns={'site_id':'sensor_name','latitude':'lat','longitude':'long'})
    else:
        print('Sensor file type not recognised: {}'.format(sensor_file_type))
        print('Should be CSV or RDATA, please amend as appropriate.')
        quit()

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


def create_daily_dataframe_for_sites(emep_ds,sens):

    # get the dates we are going to process (cutting off any dates for which we have less than 24 data points)
    dates = emep_ds.time.groupby("time.dayofyear").min()[emep_ds.time.groupby("time.dayofyear").count()==24].values     

    # get the list of sites from the sensor data
    sensor_list = sens.sensor_name

    ### create dataframe for storing model data
    #     this will start off with the format:
    #    rows: [date] [sensor name]
    # columns: [SPC Daily Mean] [SPC Daily Max]

    arrays = [dates,sensor_list]
    name_strings = ['Date','SiteID']
    dsind = pd.MultiIndex.from_product(arrays,names=name_strings)

    col_strings = ['NO2_mean','NO2_max','SO2_mean','SO2_max','NOx_mean','NOx_max','PM2.5_mean','PM2.5_max','PM10_mean','PM10_max']

    model_data = pd.DataFrame(index=dsind,columns=col_strings)

    return(model_data)


def create_full_dataframe_for_sites(emep_ds,sens):

    # get the dates we are going to process (cutting off any dates for which we have less than 24 data points)
    dates = emep_ds.time.values     

    # get the list of sites from the sensor data
    sensor_list = sens.sensor_name

    ### create dataframe for storing model data
    #     this will start off with the format:
    #    rows: [date] [sensor name]
    # columns: [SPC]

    arrays = [dates,sensor_list]
    name_strings = ['Date','SiteID']
    dsind = pd.MultiIndex.from_product(arrays,names=name_strings)

    col_strings = ['NO2','SO2','NOx','PM2.5','PM10']

    model_data = pd.DataFrame(index=dsind,columns=col_strings)

    return(model_data)


def calc_daily_max(dataframe):
    return(dataframe.groupby("time.dayofyear").max().where(dataframe.groupby("time.dayofyear").count()==24))

def calc_daily_mean(dataframe):
    return(dataframe.groupby("time.dayofyear").mean().where(dataframe.groupby("time.dayofyear").count()==24))


def calc_daily_data(model_data,NO2_data,NOX_data,SO2_data,O3_data,PM25_data,PM10_data):

    # calculate the daily max and mean values for whole domain
    max_NO2  = calc_daily_max(NO2_data)
    mean_NO2 = calc_daily_mean(NO2_data)

    max_NOX  = calc_daily_max(NOX_data)
    mean_NOX = calc_daily_mean(NOX_data)

    max_SO2  = calc_daily_max(SO2_data)
    mean_SO2 = calc_daily_mean(SO2_data)

    max_O3  = calc_daily_max(O3_data)
    mean_O3 = calc_daily_mean(O3_data)

    max_PM25  = calc_daily_max(PM25_data)
    mean_PM25 = calc_daily_mean(PM25_data)

    max_PM10  = calc_daily_max(PM10_data)
    mean_PM10 = calc_daily_mean(PM10_data)


    # loop through sites, and dates of interest, to extract the required data
    for siteID, temp_df in model_data.groupby(level=1):
        i_index = sens[sens.sensor_name==siteID].iarray.values[0]
        j_index = sens[sens.sensor_name==siteID].jarray.values[0]


        for t_index,dstamp in enumerate(temp_df.index.get_level_values(level='Date')):
            model_data.loc[(dstamp,siteID)]['NO2_mean'] = mean_NO2[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['NO2_max']  = max_NO2[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['NOx_mean'] = mean_NOX[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['NOx_max']  = max_NOX[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['SO2_mean'] = mean_SO2[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['SO2_max']  = max_SO2[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['O3_mean']  = mean_O3[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['O3_max']   = max_O3[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['PM2.5_mean'] = mean_PM25[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['PM2.5_max']  = max_PM25[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['PM10_mean'] = mean_PM10[t_index,j_index,i_index].values
            model_data.loc[(dstamp,siteID)]['PM10_max']  = max_PM10[t_index,j_index,i_index].values


    return(model_data)


def calc_full_data(model_data,NO2_data,NOX_data,SO2_data,O3_data,PM25_data,PM10_data):

    model_data = model_data.reorder_levels([1,0])
    model_data = model_data.sort_index(axis=0)

    # loop through sites, and dates of interest, to extract the required data
    for siteID  in model_data.index.levels[0]:
        i_index = sens[sens.sensor_name==siteID].iarray.values[0]
        j_index = sens[sens.sensor_name==siteID].jarray.values[0]

        model_data.loc[(siteID,),'NO2']   = NO2_data[:,j_index,i_index].values
        model_data.loc[(siteID,),'NOx']   = NOX_data[:,j_index,i_index].values
        model_data.loc[(siteID,),'SO2']   = SO2_data[:,j_index,i_index].values
        model_data.loc[(siteID,),'O3']    = O3_data[:,j_index,i_index].values
        model_data.loc[(siteID,),'PM2.5'] = PM25_data[:,j_index,i_index].values
        model_data.loc[(siteID,),'PM10']  = PM10_data[:,j_index,i_index].values


    model_data = model_data.reorder_levels([1,0])
    model_data = model_data.sort_index(axis=0)
    return(model_data)



def extract_emep_data(emep_ds,sens,daily_means):

    # create the dataframe for our site data (species list is hardcoded in this function)
    if daily_means:
        model_data = create_daily_dataframe_for_sites(emep_ds,sens)
    else:
        model_data = create_full_dataframe_for_sites(emep_ds,sens)

    # load the required data (hardcoded for the moment)
    NO2_data = emep_ds.SURF_ppb_NO2
    NOX_data = emep_ds.SURF_ppb_NO2 + emep_ds.SURF_ppb_NO
    SO2_data = emep_ds.SURF_ppb_SO2
    O3_data  = emep_ds.SURF_ppb_O3
    PM25_data = emep_ds.SURF_ug_PM25X_rh50
    PM10_data = emep_ds.SURF_ug_PM10

    if daily_means:
        model_data = calc_daily_data(model_data,NO2_data,NOX_data,SO2_data,O3_data,PM25_data,PM10_data)
    else:
        model_data = calc_full_data(model_data,NO2_data,NOX_data,SO2_data,O3_data,PM25_data,PM10_data)
    
    return(model_data)



if __name__ == '__main__':

    import argparse
    
    # read arguments from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--sensor_file","-s",help="path and name of the sensor locations file")
    parser.add_argument("--sensor_file_type","-t",help="type of sensor file: CSV or RDATA")
    parser.add_argument("--wrf_file","-w",help="path and name of the WRF output file for location information")
    parser.add_argument("--emep_file","-e",help="path and name of the EMEP output file to process")
    parser.add_argument("--out_file","-o",help="path and name of the EMEP output file to process")
    parser.add_argument("--daily_means",dest="daily_mean",action='store_true',help="flag to indicate that daily mean/max data should be output")
    parser.add_argument("--no_daily_means",dest="daily_mean",action='store_false',help="flag to indicate that raw data should be output")
    parser.set_defaults(daily_mean=True)
    

    # read arguments from the command line
    args = parser.parse_args()

    if args.sensor_file:
        sensor_file = Path(args.sensor_file)
    else:
        #sensor_file = Path("./sensors_detailed.csv")
        print('Please provide sensor file location, using --sensor_file')

    if args.sensor_file_type:
        sensor_file_type = args.sensor_file_type
    else:
        #sensor_file = Path("./sensors_detailed.csv")
        print('Please provide sensor file type, using --sensor_file_type')


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

    if not args.daily_mean:
        print('We will be outputting data for all time stamps in the original dataset')


    ### operational code

    print('loading sensor site locations')
    # open the sensor dataset, with grid locations
    sens = load_sites_and_obtain_their_grid_locations(wrf_in,sensor_file,sensor_file_type)

    print('opening EMEP file')
    # open EMEP file
    emep_ds = xr.open_dataset(emep_in)
    # indexing later will be done:
    #    emep_ds.T2C[time,jarray,iarray]  
    
    print('extracting EMEP data')
    # extract the model data for each site
    site_dataframe = extract_emep_data(emep_ds,sens,args.daily_mean)

    print('writing to csv file')
    # write the model data to csv
    site_dataframe.to_csv(out_file)

    print('finished')











