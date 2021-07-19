"""
Script for combining the individual pollen datasets into a single datafile,
while cleaning the dates and site ID's too.

@author:    Douglas Lowe and Ann Gledson
            Research IT, University of Manchester
"""


import argparse
import sys
sys.path.append("../..")
from pathlib import Path

from environmental_data_modules import PollenPostProcessor


if __name__ == '__main__':

    ### read arguments from the command line
    help_string = """Pollen datasets processing script."""
    
    parser = argparse.ArgumentParser(description=help_string)
    
    parser.add_argument("--data_dir","-d",dest="data_dir",type=str,
                        help="data directory name. Default: {}".format(PollenPostProcessor.DEFAULT_DATA_DIR))
    
    parser.add_argument("--file_template","-f",dest="file_template",type=str,
                        help="file template: Default: {}".format(PollenPostProcessor.DEFAULT_FILE_TEMPLATE))
    
    parser.add_argument("--outfile","-o",dest="out_file",type=str,
                        help="out file name: Default is None")
    parser.add_argument("--outfile_suffix","-s",dest="out_suffix",type=str,
                        help="out file suffix, to use with file template: Default is None")

    parser.add_argument("--pollen_list","-p",dest="pollenlist",type=list,
                        help="list of pollen species to process: Default: [{}]".format(','.join(PollenPostProcessor.DEFAULT_POLLEN)))

    # Log verbose-ness
    parser.add_argument("--verbose", "-v", type=int,
                        help="Level of output for debugging (Default: {} (0 = no verbose output))".format(
                            str(PollenPostProcessor.DEFAULT_VERBOSE)))

    ### Process inputs
    args = parser.parse_args()

    if args.data_dir:
        data_dir = args.data_dir
    else:
        data_dir = PollenPostProcessor.DEFAULT_DATA_DIR
    
    if args.file_template:
        file_template = args.file_template
    else:
        file_template = PollenPostProcessor.DEFAULT_FILE_TEMPLATE
    
    if args.out_file:
        out_file = args.out_file
    else:
        out_file = None
    
    if args.out_suffix:
        out_suffix = args.out_suffix
    else:
        out_suffix = None
        
    if args.pollenlist:
        pollenlist = args.pollenlist
    else:
        pollenlist = PollenPostProcessor.DEFAULT_POLLEN 

    
    if args.verbose:
        verbose = max(args.verbose, 0)
        print('output verbose level: {}'.format(verbose))
    else:
        print('No verbose flag provided, so using default: {}'.format(str(PollenPostProcessor.DEFAULT_VERBOSE)))
        verbose = PollenPostProcessor.DEFAULT_VERBOSE


    pollen_processor = PollenPostProcessor(verbose=verbose)
    
    pollendata = pollen_processor.process(file_template=file_template,outfile=out_file,outsuffix=out_suffix,
                             data_directory=data_dir,pollen_species=pollenlist,save_to_csv=True)



