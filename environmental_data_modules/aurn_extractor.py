from abc import ABCMeta, abstractmethod

from .environment_module import EnvironmentModule


class AurnExtractor(EnvironmentModule):
    __metaclass__ = ABCMeta


    def __init__(self, out_dir=EnvironmentModule.DEFAULT_OUT_DIR, verbose=EnvironmentModule.DEFAULT_VERBOSE):
        super(AurnExtractor, self).__init__(out_dir, verbose)
        self._headstring = None
        self._cols_specific = []
        self._head_string = 'Rain gauge daily data{} for date range: {} to {}\n'
        self._cols_specific = []
        self._file_out = '{}/rain{}{}.csv'.format(self.out_dir, '{}', '{}')


    def _get_settings(self, outfile_suffix=EnvironmentModule.DEFAULT_OUT_FILE_SUFFIX):
        filename = self._file_out.format('', outfile_suffix)
        headstring = self._head_string.format('', self.date_range[0], self.date_range[1])

        return {
            'fname': filename,
            'headstring': headstring,
            'columnstring': ','.join(self.get_all_cols()) + '\n'
        }



    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        settings = self._get_settings(outfile_suffix)
        return self._extract_data(settings, self._save_to_file)


    def _extract_data(self, settings, save_to_file=True):
        if self.verbose >= 1:
            print('extracting data for {}')
        if self.verbose > 1:
            print('using extraction dict: {}')
            print('using settings: {}'.format(settings))
        datadata = self._perform_extraction()

        return datadata

    def _perform_extraction(self):
        datadata = {}

        return datadata

    def _save_to_file(self, data_result, settings):
        print('saving to file: {}'.format(settings['fname']))

        with open(settings['fname'], 'w') as dfile:
            dfile.write(settings['headstring'])
            dfile.write(settings['columnstring'])

            for data in data_result.values():

                dfile.write('\n')
        if self.verbose > 0:
            print('saved')
