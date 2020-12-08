from abc import ABCMeta

class MetModule(object):
    """
        Abstract Class used for extracting and post-processing data that has been extracted from MEDMI server.

    """
    __metaclass__ = ABCMeta

    COLUMNS_BASE = ['timestamp', 'site_id']

    def __init__(self):
        """ Initialise instance of the MetModule class.
            Initialises the private class variables with hard-coded / default values

            Returns:
                Initialised instance of subclass of MetModule

        """
        self._columns_base = MetModule.COLUMNS_BASE
        self._columns_specific = []

    def get_all_column_headers(self):
        """ Get all the column headers for the output CSV file

            Returns:
                List of strings: columns headers used for Met extraction/post-processing output files

        """
        return self._columns_base + self._columns_specific
