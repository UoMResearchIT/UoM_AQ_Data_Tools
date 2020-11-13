from abc import ABCMeta, abstractmethod

class MetModule(object):
    __metaclass__ = ABCMeta
    """
        Abstract Class used for extracting and post-processing data that has been extracted from MEDMI server.
    """
    COLUMNS_BASE = ['date', 'siteID']

    def __init__(self):
        self._columns_base = MetModule.COLUMNS_BASE
        self._columns_specific = []

    def get_all_columns(self):
        return self._columns_base + self._columns_specific
