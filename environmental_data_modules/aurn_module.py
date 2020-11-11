try:
    import wget
    import pyreadr
    from pathlib import Path
    from abc import ABCMeta, abstractmethod
except:
    pass


class AurnModule(object):
    INDEX_EXTRACTED = 'index'
    SITE_ID_EXTRACTED = 'siteID'
    SITE_ID_NEW = 'SiteID'
    SITE_ID_AURN_METADATA = 'site_id'
    EXTRACTED_FILE_COLS = ['Date', 'SiteID', 'O3', 'NO2', 'NOXasNO2', 'SO2', 'PM2.5', 'PM10']

    # Define defaults
    DEFAULT_METADATA_FILE = "AURN_metadata.RData"
    DEFAULT_METADATA_URL = 'https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData'
    DEFAULT_SITE_LIST = None

    def __init__(self, metadata_filename=DEFAULT_METADATA_FILE, metadata_url=DEFAULT_METADATA_URL):
        self._metadata = self.load_metadata(metadata_filename, metadata_url)
        self._site_list = AurnModule.DEFAULT_SITE_LIST

    @property
    def metadata(self):
        return self._metadata

    @property
    def site_list(self):
        return self._site_list

    @site_list.setter
    def site_list(self, site_list):
        all_sites = self.metadata['AURN_metadata'][AurnModule.SITE_ID_AURN_METADATA].unique()
        # get list of sites to process extract_site_data
        if site_list is None:
            self._site_list = all_sites
        else:
            try:
                site_list = set(list(site_list))
            except Exception:
                raise ValueError('Site list must be a list. Input: {}'.format(site_list))
            error_sites = set(site_list) - set(all_sites)
            assert len(error_sites) == 0, \
                "Each site must be contained in available sites (from Aurn metadata: {}. Error sites: {}".format(
                    all_sites, str(error_sites))
            self._site_list = site_list


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
