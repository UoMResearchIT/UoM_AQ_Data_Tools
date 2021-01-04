from abc import ABCMeta, abstractmethod
try:
    import wget
    import pyreadr
    from pathlib import Path
except:
    pass


class AurnModule(object):
    """
        Abstract Class used for extracting and post-processing data that has been extracted from AURN server.

    """
    __metaclass__ = ABCMeta

    SPECIES_LIST_EXTRACTED = ['O3', 'PM10', 'PM2.5', 'NO2', 'NOXasNO2', 'SO2']
    INDEX_EXTRACTED = 'index'
#    SITE_ID_EXTRACTED = 'siteID'
    SITE_ID_AURN_METADATA = 'site_id'
    DATE_EXTRACTED = 'date'

#    SITE_ID_NEW = 'SiteID'
#    DATE_NEW = 'Date'
    TIMESTAMP_STRING = 'timestamp'
    SITE_STRING = 'site_id'
    NEW_FILE_COLS = [TIMESTAMP_STRING, SITE_STRING] + SPECIES_LIST_EXTRACTED

    # Define defaults
    DEFAULT_METADATA_FILE = "AURN_metadata.RData"
    DEFAULT_DOWNLOAD_RDATA_URL = "https://uk-air.defra.gov.uk/openair/R_data/"
    DEFAULT_METADATA_URL = '{}/{}'.format(DEFAULT_DOWNLOAD_RDATA_URL, DEFAULT_METADATA_FILE)
    DEFAULT_SITE_LIST = None

    def __init__(self, metadata_filename=DEFAULT_METADATA_FILE, metadata_url=DEFAULT_METADATA_URL):
        """ Initialise instance of the MetModule class.
            Initialises the private class variables with hard-coded / default values

            Args:
                metadata_filename: filename of the metadata used in Aurn data extraction
                metadata_url: alternative source of AURN metadata, if metadata_filename is None

            Returns:
                Initialised instance of subclass of MetModule

        """
        self._timestamp_string = AurnModule.TIMESTAMP_STRING
        self._site_string = AurnModule.SITE_STRING
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
        """ Set the site_list property.
            If site_list is None, will use the site lists from the AURN metadata property
            Otherwise, will check that each site ID in site_list is in the list of all site IDs in the AURN Metadata

            Dependencies:
                self.metadata

            Args:
                site_list: list of site identifiers (strings or numbers)

            Returns:
                None

        """
        all_sites = self.metadata['AURN_metadata'][AurnModule.SITE_ID_AURN_METADATA].unique()
        # get list of sites to process extract_site_data
        if site_list is None:
            self._site_list = all_sites
        else:
            try:
                site_list = set(list(site_list))
            except Exception:
                raise TypeError('Site list must be a list. Input: {}'.format(site_list))
            error_sites = set(site_list) - set(all_sites)
            assert len(error_sites) == 0, ValueError(
                "Each site must be contained in available sites (from Aurn metadata: {}. Error sites: {}".format(
                    all_sites, str(error_sites)))
            self._site_list = site_list


    def load_metadata(self, filename=None, alt_url=None):
        """ Load the AURN metadata from file or URL.
            If filename is None, will use the metadata stored at the URL: alt_url
            Otherwise, will load from URL: alt_url

            Dependencies:
                Path

            Args:
                filename: (string) Valid file name of existing AURN metadata R file, or None
                alt_url: (string) Valid URL pointing to AURN metadata downloadable source, or None

            Returns:
                None

        """
        # Has a filname been entered and does the file exist?
        if filename is not None and Path(filename).is_file():
            print("Metadata file {} exists so will use this".format(filename))
            filename = Path(filename)
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
            raise ValueError('Invalid metadata filename and no url alternative provided')

        # Read the RData file into a Pandas dataframe
        try:
            print('Reading filename {} into dataframe.'.format(filename.name))
            return pyreadr.read_r(filename.name)
        except Exception as err:
            raise ValueError('Error reading into dataframe from R file: {} . {}'.format(filename, err))
