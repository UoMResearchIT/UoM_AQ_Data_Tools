#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 16:22:25 2020

Script for determining the AURN site addresses, and extracting postcode information,
using the googlemaps API.

This requires the following csv file as input:
    aurn_measurement_sites.csv
And creates this csv file:
    aurn_measurement_sites_addresses_postcodes.csv

To run the script use:
    python aurn_address_postcode_calculator.py --key [Google API key]


@author: mbessdl2
"""

import json
import requests
import pandas as pd
import re


#%% function for obtaining the postcode area from google api

def find_google_address(longitude,latitude,key):

    URL = 'https://maps.googleapis.com/maps/api/geocode/json?'
    
    PARAMS = {'latlng':str(latitude)+','+str(longitude),'key':key}

    rtest = requests.get(url = URL, params = PARAMS)

    postcode_dict = json.loads(rtest.content.decode("utf-8"))

    if(postcode_dict['results']):
        address = 'null'
        for results in postcode_dict['results']:
            if(results['formatted_address']):
                address = results['formatted_address']
                postcode_full = re.findall('\w{1,2}\d{1,2}\s+\d{1,2}\w{1,2}',address)
                postcode_partial = re.findall('\w{1,2}\d{1,2}',address)
                if(len(postcode_full)!=0 or len(postcode_partial)!=0):
                    break
    else:
        address = 'null'

    return(address)


def extract_postcode(address):
        
    postcode = 'error'
    
    if(address!='null'):
        try:
            postcode_full = re.findall('[A-Za-z]{1,2}\d{1,2}\s+\d{1,2}[A-Za-z]{1,2}',address)[0]
            postcode = re.subn('\d*\s+\d*[A-Za-z]*','',postcode_full)[0]
        except IndexError:
            try:
                postcode_full = re.findall('[A-Za-z]{1,2}\d{1,2}[A-Za-z]{1}\s+\d{1,2}[A-Za-z]{1,2}',address)[0]
                postcode = re.subn('\d*[A-Za-z]*\s+\d*[A-Za-z]*','',postcode_full)[0]
            except IndexError:                
                try:
                    postcode_full = re.findall('[A-Za-z]{1,2}\d{1,2}',address)[0]
                    postcode = re.subn('\d*','',postcode_full)[0]
                except:
                    postcode = 'null'
            except:
                postcode = 'null'
        except:
            postcode = 'null'
    else:
        postcode = 'null'

    return(postcode)


if __name__ == '__main__':

    import argparse
    
    # initiate the parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--key", "-k", help="google API key")
    
    # read arguments from the command line
    args = parser.parse_args()

    if args.key:
        google_key = args.key
    else:
        print('please pass the google API key, using --key')
        
    
    print('reading station data')
    site_metadata = pd.read_csv('aurn_measurement_sites.csv')
    site_metadata = site_metadata.drop(['EU Site ID', 'EMEP Site ID',
           'Zone', 'Start Date', 'End Date',
           'Northing', 'Easting', 'Altitude (m)', 'Networks',
           'AURN Pollutants Measured', 'Site Description'],axis=1)
    
    print('obtaining addresses from google API (will take some time)')
    site_metadata['Address'] = site_metadata.apply(lambda row: find_google_address(row.Longitude,row.Latitude,google_key), axis=1)
    print('extracting postcode data')
    site_metadata['Postcode'] = site_metadata.apply(lambda row: extract_postcode(row.Address), axis=1)
    
    print('writing output file')
    site_metadata.to_csv('aurn_measurement_sites_addresses_postcodes.csv',na_rep='NaN')

