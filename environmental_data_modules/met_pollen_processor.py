import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path

from environmental_data_modules import PostProcessor, MetModule, DateRangeProcessor


class PollenPostProcessor(PostProcessor, MetModule, DateRangeProcessor):
    """
        Class used for post-processing the pollen data extracted from the MEDMI server.
    """

    # Define 'absolute' constants
    DATE_CALCS_FORMAT = '%Y-%m-%d %H:%M:%S'
    POLLEN_SITE_STRING = ' [POLLEN]'

   # Define default constants
    DEFAULT_DATA_DIR = 'full_data'
    DEFAULT_FILE_TEMPLATE = 'Met_extracted_pollen-{}_full.csv'
    DEFAULT_POLLEN = ['alnus','ambrosia','artemisia','betula','corylus', \
              'fraxinus','platanus','poaceae','quercus','salix', \
              'ulmus','urtica']

    def __init__(self, out_dir='./', verbose=PostProcessor.DEFAULT_VERBOSE):
        """ Initialise instance of the PollenPostProcessor class.
        
        """
        super(PollenPostProcessor, self).__init__(out_dir, verbose)
        MetModule.__init__(self)
        DateRangeProcessor.__init__(self)


    def process(self, file_template=DEFAULT_FILE_TEMPLATE, 
                outfile=None, outsuffix=None, data_directory=DEFAULT_DATA_DIR,
                pollen_species=DEFAULT_POLLEN, save_to_csv=PostProcessor.DEFAULT_SAVE_TO_CSV):
        """ Post process for the pollen data from MEDMI.
        
        This is simpler than the Met or Pollution post processor.
        We tidy up the date strings, extend the site ID's to include
        the [POLLEN] string, and combine the individual pollen datasets
        into a single dataset.
        
        Args:
            in_files:
            pollen_species:
            save_to_csv: 
            outfile_suffix:
        
        
        Returns:
            pollen_data_daily: daily dataset, for all pollen species and sites, as pandas.Dataframe
                Required MultiIndex:
                        'time_stamp'  (datetime object): date (only) (e.g. 2017-06-01)
                        'sensor_name'          (string): ID string for site (e.g. '3 [POLLEN]')
                Required Columns:
                        '[pollen]'              (float): pollen count for given pollen species                
        
        """
               
        work_dir = Path(data_directory)
        assert work_dir.is_dir()
        if outfile:
            outfilepath = work_dir.joinpath(outfile)
        elif outsuffix:
            outfilepath = work_dir.joinpath(file_template.format(outsuffix))
        else:
            print('Must specify outfile or outsuffix')
            return
        
        
        index_array = [self._timestamp_string,self._site_string]
        
        for pollen in pollen_species:
            work_file = work_dir.joinpath(file_template.format(pollen))
            if work_file.is_file():
                temp_dataset = pd.read_csv(work_file,skiprows=1)
                temp_dataset[self._timestamp_string] = temp_dataset[self._timestamp_string].apply(self.parse_calcs_date_only)
                temp_dataset[self._site_string] = temp_dataset[self._site_string].apply(self.rename_siteid)                
                temp_dataset = temp_dataset.set_index(index_array)
                
                try:
                    pollen_dataset = pd.merge(pollen_dataset,temp_dataset,how='outer',on=index_array)
                except:
                    pollen_dataset = temp_dataset.copy(deep=True)
                
        
        if save_to_csv:
            pollen_dataset.to_csv(outfilepath, index=True, header=True)
    
        return pollen_dataset
    
    ### function for creating date objects
    def parse_calcs_date_only(self, date_in):
        tempdate = datetime.strptime(date_in, self.DATE_CALCS_FORMAT)
        return tempdate.date()

    def rename_siteid(self, data_in):
        return str(data_in)+self.POLLEN_SITE_STRING

    
