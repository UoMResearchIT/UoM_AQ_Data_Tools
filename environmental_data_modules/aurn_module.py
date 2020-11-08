import wget
import pyreadr
from pathlib import Path
from abc import ABCMeta, abstractmethod


class AurnModule(object):
    __metaclass__ = ABCMeta
    EXTRACTED_FILE_INDEX = 'index'
    EXTRACTED_FILE_COLS = ['Date', 'SiteID', 'O3', 'NO2', 'NOXasNO2', 'SO2', 'PM2.5','PM10']

    # Define defaults
    DEFAULT_METADATA_FILE = "AURN_metadata.RData"
    DEFAULT_METADATA_URL = 'https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData'

    def __init__(self, metadata_filename=DEFAULT_METADATA_FILE, metadata_url=DEFAULT_METADATA_URL):
        self._metadata = self.load_metadata(metadata_filename, metadata_url)

    @property
    def metadata(self):
        return self._metadata

    def load_metadata(self, filename, alt_url=None):
        # Does the file exist?
        filename = Path(filename)
        if filename.is_file():
            print("Metadata file {} already exists so will use this".format(filename))
        elif alt_url:
            # Does the URL alternative exist and does it work
            print("Downloading data file using url {}".format(alt_url))
            try:
                filename = Path(wget.download(alt_url))
                print('\nMetadata file loaded from url')
            except Exception as err:
                raise ValueError('Error obtaining metadata file from {}. {}'.format(alt_url, err))
        else:
            # Neither works
            raise ValueError('Metadata filename does not exist and no url alternative provided')

        # Read the RData file into a Pandas dataframe
        try:
            print('Reading filename {} into dataframe.'.format(filename.name))
            return pyreadr.read_r(filename.name)
        except Exception as err:
            raise ValueError('Error reading into dataframe from R file: {} . {}'.format(filename, err))
