from abc import ABCMeta, abstractmethod
import os, errno
from datetime import datetime


class EnvironmentModule:
    """
        Abstract class, parent of classes used for extracting data from AURN/MEDMI servers and post-processing.
    """
    __metaclass__ = ABCMeta

    # Define 'absolute' constants
    UK_LATITUDES = [48., 60.]
    UK_LONGITUDES = [-11., 3.]
    LONGITUDE_RANGE = [-180., 360.]
    LATITUDE_RANGE = [-90., 90.]

    # Define defaults
    DEFAULT_OUT_FILE_SUFFIX = ''
    DEFAULT_OUT_DIR = 'out_dir'
    DEFAULT_START_DATE = datetime(2016, 1, 1, 0)
    DEFAULT_END_DATE = datetime(2019, 12, 31, 23)
    DEFAULT_DATE_RANGE = [DEFAULT_START_DATE, DEFAULT_END_DATE]
    DEFAULT_VERBOSE = 0


    def __init__(self, dir_name=DEFAULT_OUT_DIR, verbose=DEFAULT_VERBOSE):
        # Todo - What are the verbose limits? 0..n
        """ Initialise instance of the EnvironmentModule class.
            Args:
                dir_name: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of subclass of EnvironmentModule

        """
        self.out_dir = dir_name
        self.verbose = verbose
        self.__date_range = EnvironmentModule.DEFAULT_DATE_RANGE
        self._base_file_out = None
        self._file_out = None
        self._outfile_suffix = EnvironmentModule.DEFAULT_OUT_FILE_SUFFIX



    @staticmethod
    def all_subclasses(cls):
        return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in EnvironmentModule.all_subclasses(c)])

    @property
    def outfile_suffix_string(self):
        return '_' + self._outfile_suffix if len(self._outfile_suffix.strip()) > 0 else ''

    @property
    def base_file_out(self):
        return self._base_file_out

    @property
    def date_range(self):
        return self.__date_range

    @date_range.setter
    def date_range(self, range):
        # Example input: ['2017-1-1_0', '2017-06-30_23']
        if isinstance(range[0], datetime):
            datetime_1 = range[0]
        else:
            raise ValueError('Start date is not in datetime format')
        if isinstance(range[1], datetime):
            datetime_2 = range[1]
        else:
            raise ValueError('End date is not in datetime format')
        if datetime_1 >= datetime_2:
            raise ValueError('Start date is not earlier than end date.')
        self.__date_range = [datetime_1, datetime_2]
        

    @property
    def out_dir(self):
        return self.__out_dir

    @out_dir.setter
    def out_dir(self, dir_name):
        try:
            dir_name = str(dir_name)
        except ValueError as err:
            raise err
        try:
            os.makedirs(dir_name)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise ValueError("Directory name {} cannot be created.".format(dir_name))
        self.__out_dir = dir_name

    @property
    def verbose(self):
        return self.__verbose

    @verbose.setter
    def verbose(self, verbose):
        try:
            verbose = int(verbose)
        except ValueError as err:
            raise err
        self.__verbose = max(0, verbose)

    @property
    def file_out(self):
        return self._file_out

    @file_out.setter
    def file_out(self, filename):
        try:
            filename = str(filename)
        except ValueError as err:
            raise err
        #Todo check if filename is createable
        self._file_out = filename

