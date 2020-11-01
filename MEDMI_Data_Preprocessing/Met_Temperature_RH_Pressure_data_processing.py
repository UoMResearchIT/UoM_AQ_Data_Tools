import argparse
import sys
sys.path.append("..")
from environmental_data_workflows import MetPostProcessor


if __name__ == '__main__':
    # Todo add arguments/processing

    impute_data = True
    print_stats = True
    #file_in = 'inputs/temp_extras-rel_hum-pressure-dewpoint_June2017.csv'
    file_in = 'inputs/temp_extras-rel_hum-pressure-dewpoint_midlands.csv'
    file_out = 'daily_mean_max_temp_RH_pres.csv'
    station_drop_list = [117]
    min_temperature = -20
    date_range = ['2017-01-01 00', '2019-06-30 23']
    reference_station_number = 5
    verbose = 2
    min_years = 0.04 #1
    reference_num_years = 0.0625 #3.5
    out_dir = 'met_postprocessing'

    post_processor = MetPostProcessor(out_dir, verbose)
    post_processor.post_process(file_in, date_range=date_range, file_out=file_out, station_drop_list=station_drop_list,
                                min_temperature=min_temperature, reference_station_number=reference_station_number,
                                min_years=min_years, reference_num_years=reference_num_years,
                                impute_data=impute_data, print_stats=print_stats)
