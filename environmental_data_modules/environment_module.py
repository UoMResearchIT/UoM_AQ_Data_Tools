from abc import ABCMeta, abstractmethod
import os, errno

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

