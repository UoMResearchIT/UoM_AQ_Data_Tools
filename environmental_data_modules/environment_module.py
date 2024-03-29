from abc import ABCMeta, abstractmethod
import os, errno

class EnvironmentModule:
    """
        Abstract class, parent of classes used for extracting data from AURN/MEDMI servers and post-processing.
    """
    __metaclass__ = ABCMeta

    # Define defaults
    DEFAULT_OUT_FILE_SUFFIX = ''
    DEFAULT_OUT_DIR = 'out_dir'
    DEFAULT_VERBOSE = 0


    def __init__(self, out_dir=DEFAULT_OUT_DIR, verbose=DEFAULT_VERBOSE):
        """ Initialise instance of the EnvironmentModule class.
            Args:
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output. Zero = no output

            Returns:
                Initialised instance of subclass of EnvironmentModule

        """
        self.out_dir = out_dir
        self.verbose = verbose
        self._base_file_out = None
        self._file_out = None
        self._outfile_suffix = EnvironmentModule.DEFAULT_OUT_FILE_SUFFIX

    ### Static methods ###

    @staticmethod
    def all_subclasses(cls):
        return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in EnvironmentModule.all_subclasses(c)])

    ### Class properties: Get/Sets ###

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
        assert isinstance(dir_name, str), ValueError("dir_name is not a valid string")
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
        assert isinstance(verbose, int), ValueError("verbose is not a valid integer")
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

