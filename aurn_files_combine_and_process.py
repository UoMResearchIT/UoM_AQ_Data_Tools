#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:22:12 2020

Script for reading in AURN measurement data, and outputting a combined, formatted, data file.

It is run using the command:
	python aurn_files_combine_and_process.py -i [input directory] -o [output directory] -s [AURN station file]


AURN data taken from this website:
    https://uk-air.defra.gov.uk/data/data_selector
    (select 'Search Hourly Networks', then 'Daily Mean' or 'Daily Max', and use AURN datasets)

Inputs:
    files of naming format [string].csv

    Example 1st 3 lines of headers for inputs:
        Daily Mean data	    supplied by UK-air on 28/01/2020
        All Data GMT hour ending
        Status: V=Verified P=Provisionaly Verified N=Not Verified S=Suspect

    Example lines 4 & 5 of headers for inputs:
        ,"Manchester Piccadilly",,"Middlesbrough",
        Date,"Ozone",Status,"Ozone",Status

    Lines 4 & 5 will be used as headers for dataframe, with a 3rd line (6) added
    to indicate 'mean' or 'max' (following information in 1st line of header)

    NOTE:
    Some input files downloaded from AURN seem to be missing information from the 1st part
    of the last column. This script will fix this problem when creating the intermediate
    data files.
    

Intermediate files:
    These will have the gaps in lines 4 & 5 filled in, and 'Status' header will be given
    the species information. Line 6 of the header will be added. Status flag will be 
    a simple V or N. If any lines are missing the last column this will be added (as an empty
    column).
    
    File names will be: [string].csv.headers.csv

Station information input file:
	This contains information on the location and postcode for each of the AURN stations.
	The default name assumed is './aurn_measurement_sites_addresses_postcodes.csv', 
	but another filename can be specified, using the '-s' flag, as described above.

Outputs:
	A csv file containing the combined, sorted, data. This will contain the columns:
	date, location, [datasets], Address, PostCode, lati, loni
	


@author: mbessdl2
"""



import numpy as np
import pandas as pd
import glob
import re
import os


def delete_old_temp_files(input_dir):
    removelist = glob.glob(input_dir+"/*.csv.headers.csv")

    for rfile in removelist:
        os.remove(rfile)


# functions for making header strings more readable
def shorten_chem_spec(spc_string):
    spc_dict = {}
    spc_dict['Nitrogen oxides as nitrogen dioxide'] = 'NOx'
    spc_dict['Nitrogen dioxide'] = 'NO2'
    spc_dict['Ozone'] = 'O3'
    spc_dict['PM10 particulate matter (Hourly measured)'] = 'PM10'
    spc_dict['PM2.5 particulate matter (Hourly measured)'] = 'PM2.5'
    spc_dict['Sulphur dioxide'] = 'SO2'

    spc_string = spc_string.strip('\"')

    outstring = [value for key, value in spc_dict.items() if spc_string in key]

    return(outstring)


# function for correcting the headers & column count of the input files.
#    this creates a new set of (temporary) input files, with the endings: *.headers.csv
def correct_headers_and_column_count(input_dir):
    filelist = glob.glob(input_dir+"/*.csv")
    
    end_string = re.compile('end\s*\n',re.IGNORECASE)
    mean_string = re.compile('mean',re.IGNORECASE)
    max_string  = re.compile('maximum',re.IGNORECASE)
    missing_last_column_string = re.compile('No data\n',re.IGNORECASE)
    
    status_string = re.compile(' \(FIDAS\)| \(TEOM FDMS\)| \(BAM\)| \(Ref\.eq\)| ugm-3')
    
    
    for filename in filelist:
        
        # skip the output file, if it exists
        if(re.match('AURN_AQ_data',filename)):
            continue
        
        with open(filename,'r') as f:
            filelines = f.readlines()
    
        s=","
    
        ### deal with the missing spaces for the location list
        fileareas = filelines[3].rstrip().split(',')    
        
        fileareas[0] = 'Location'   # create the index value for this row
        
        for ind, area in enumerate(fileareas):
            if(area == ''):
                fileareas[ind] = fileareas[ind-1]
    
        filelines[3] = s.join(fileareas)+'\n'
        
        
        
        ### create header for labelling data type (values or status)
        ### also create a header list with the species name
        fileinfo = filelines[4].rstrip().split(',')
        filespc  = filelines[4].rstrip().split(',')
        
        for ind, info in enumerate(fileinfo):
            if(info == 'Status'):
                filespc[ind]  = filespc[ind-1]
            else:
                fileinfo[ind] = 'Value'
        
        fileinfo[0] = 'Type'
        
        filelines[4] = s.join(fileinfo)+'\n'
    
    
        ### rename the species (shorten them), and add mean / max label info
        filespc = [shorten_chem_spec(x) for x in filespc]
            
        
        if(mean_string.findall(filelines[0])):
            head_string = '_mean'
        elif(max_string.findall(filelines[0])):
            head_string = '_max'
        else:
            head_string = '_error'
        
        for ind, label in enumerate(filespc):
            if(len(label)>0):
                filespc[ind] = label[0]+head_string
        filespc[0] = 'Species'
        
        filelines[3] = filelines[3] + s.join(filespc)+'\n'
    
        
        ### delete the last line (END)
        for ind, line in enumerate(filelines):
            if(end_string.match(line)):
                filelines[ind]='\n'
        
        ### convert the status strings to simple V/P/N/S
        for ind, line in enumerate(filelines):
            filelines[ind]=status_string.subn('',line)[0]
        
        ### deal with any lines which are missing the last column of data
        for ind, line in enumerate(filelines):
            filelines[ind]=missing_last_column_string.subn('No data,N\n',line)[0]
        
        ### write the data back to a temporary file
        with open(filename+'.headers.csv', 'w') as f:
            for item in filelines:
                f.write("%s" % item)
        
    

# function for reading in all data from the temporary files
def read_aurn_data(input_dir):
    filelist = glob.glob(input_dir+"/*.headers.csv")
    
    data_all = pd.DataFrame()
    
    data_dict = {}
    
    for filename in filelist:
        
        
        print(filename)
        data_dict[filename] = pd.read_csv(filename,header=[3,4,5],sep=',',
                                         index_col=0,na_values='No data')
        
        data_all = pd.concat([data_all,data_dict[filename]],axis=1)
    
    
    data_all.columns.names = ('location','variable','stat')
    data_all.index.name = ('date')

    return(data_all)


# function for breaking up & formatting of the AURN data
def extract_aurn_values_data(data_all):
    # this breaks up the dataset by 'stat' (value and status)
    for stat, new_df in data_all.groupby(level=2,axis=1):
        if(stat=='Status'):
            status_data = new_df
        elif(stat=='Value'):
            value_data = new_df
    
    # drop the unneeded 'stat' index (so that the boolean mask will work)
    status_data.columns = status_data.columns.droplevel(level=2)
    value_data.columns = value_data.columns.droplevel(level=2)
    
    # this sets any values which haven't been validated to NaN
    bool_mask = status_data != 'V'
    value_data[bool_mask] = np.nan
    
    # drop all dates without any data (so, non-validated periods)
    value_data = value_data.dropna(how='all')
        
    # first unstack puts all columns into the index, then second changes changes
    #   the variables (level=1) back into columns
    value_data = value_data.unstack().unstack(level=1)
    # swaps the two remaining indexes
    value_data = value_data.swaplevel()
    # sort the data by the first (date) index
    value_data = value_data.sort_index(axis=0,level=0)

    return(value_data)
    

# function for loading the AURN site address and postcode information
def add_address_data(value_data,aurn_file):
    
    #aurn_file = 'aurn_measurement_sites_addresses_postcodes.csv'
    
    aurn_postcodes = pd.read_csv(aurn_file,index_col=0)
    aurn_postcodes = aurn_postcodes.set_index('Site Name')
    
    # create new columns
    value_data['Address']=''
    value_data['PostCode']=''
    value_data['lati']=np.nan
    value_data['loni']=np.nan
    
    # reverse the order of 'date' 'location' in index
    value_data = value_data.swaplevel()
    
    for station in value_data.index.levels[0]:
        value_data.loc[(station,),'Address']  = aurn_postcodes.loc[station,'Address']
        value_data.loc[(station,),'PostCode'] = aurn_postcodes.loc[station,'Postcode']
        value_data.loc[(station,),'lati']     = aurn_postcodes.loc[station,'Latitude']
        value_data.loc[(station,),'loni']     = aurn_postcodes.loc[station,'Longitude']
    

    # restore the order of 'date' 'location' in index
    value_data = value_data.swaplevel()

    
    return(value_data)


# function for writing out the values dataframe
def write_out_data(value_data,output_dir):

    first_date = value_data.index.get_level_values('date').min()
    last_date = value_data.index.get_level_values('date').max()

    value_data.to_csv(output_dir+'/AURN_AQ_data_'+first_date+'_'+last_date+'.csv',na_rep='NaN')





if __name__ == '__main__':

    import argparse
    
    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", "-i", help="input directory, containing AURN .csv files")
    parser.add_argument("--output_dir", "-o", help="output directory")
    parser.add_argument("--station_file", "-s", help="AURN station information file")
    
    # read arguments from the command line
    args = parser.parse_args()

    if args.input_dir:
        input_dir = args.input_dir
    else:
        print('setting input directory to ./')
        input_dir = './'
    
    if args.output_dir:
        output_dir = args.output_dir
    else:
        print('setting output directory to ./')
        output_dir = './'
    
    if args.station_file:
        station_file = args.station_file
    else:
        print('assuming AURN station file is ./aurn_measurement_sites_addresses_postcodes.csv')
        station_file = './aurn_measurement_sites_addresses_postcodes.csv'
    
    print('deleting any old temporary csv files')
    delete_old_temp_files(input_dir)
    print('creating new temporary csv files')
    correct_headers_and_column_count(input_dir)
    
    print('loading csv files into dataframe')
    aurn_data = read_aurn_data(input_dir)
    print('pulling out values, dropping non-validated values')
    aurn_data = extract_aurn_values_data(aurn_data)
    print('adding the station address information')
    aurn_data = add_address_data(aurn_data,station_file)
    print('writing out dataset')
    write_out_data(aurn_data,output_dir)






