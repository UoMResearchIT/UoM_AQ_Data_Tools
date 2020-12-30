try:
    import os, errno
    import pandas as pd
    import numpy as np
    from pathlib import Path

    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    from scipy import stats

    sns.set_theme()
except:
    pass

daily_file = 'aurn_daily_correlation_stats.csv'
hourly_file = 'aurn_hourly_correlation_stats.csv'

path_head = 'aurn_stats/test_stats_{}_{}'

chemical_scenarios = ['NO2','NOXasNO2','O3','PM10','PM2.5','SO2','combined','combined_and_so2']
chemical_dictionary = {'NO2':'individual','NOXasNO2':'individual','O3':'individual','PM10':'individual',
                       'PM2.5':'individual','SO2':'individual','combined':'combined','combined_and_so2':'combined_and_so2'}
chemical_mixes = ['individual','combined','combined_and_so2']

sample_scenarios = ['random','random_with_EMEP','startloss','startloss_with_EMEP']

stats_list = ['kendalltau_corr','pearsonr_corr','slope']
index_list = ['site_id','spc']
col_hourly_list = stats_list + index_list
col_daily_list = col_hourly_list + ['stat']


chemical_pairs = {'NO2':'NOx','NOXasNO2':'NOx','O3':'SO2 O3','PM10':'PM','PM2.5':'PM','SO2':'SO2 O3'}





def read_datafile(full_data,filepath,index_list_internal,col_list_internal,scenario_string,chemical_string,chem_pairs_internal):
    try:
        temp_data = pd.read_csv(filepath,usecols=col_list_internal)
        temp_data['data loss scenario'] = scenario_string
        temp_data['chemical mix'] = chemical_string
        temp_data['chemical group'] = [chem_pairs_internal[chemkey] for chemkey in temp_data['spc']]
        temp_data = temp_data.rename(columns={'spc':'chemical species'})
        return(full_data.append(temp_data))
    except:
        print('could not process file {}'.format(filepath))
        return(full_data)


daily_data  = pd.DataFrame()
hourly_data = pd.DataFrame()


for scenario in sample_scenarios:
    for chemical in chemical_scenarios:
        base_dir = Path(path_head.format(chemical,scenario))

        if base_dir.is_dir():
            hourly_path = base_dir.joinpath(hourly_file) 
            daily_path  = base_dir.joinpath(daily_file)
            daily_data = read_datafile(daily_data,daily_path,index_list,col_daily_list,scenario,chemical_dictionary[chemical],chemical_pairs)
            hourly_data = read_datafile(hourly_data,hourly_path,index_list,col_hourly_list,scenario,chemical_dictionary[chemical],chemical_pairs)
        else:
            print('{} directory does not exist'.format(base_dir))

daily_mean_data = daily_data[daily_data['stat']=='mean']
daily_max_data = daily_data[daily_data['stat']=='max']
daily_mean_data.set_index(['data loss scenario','chemical mix','chemical species'],inplace=True)
daily_max_data.set_index(['data loss scenario','chemical mix','chemical species'],inplace=True)
hourly_data.set_index(['data loss scenario','chemical mix','chemical species'],inplace=True)

daily_mean_data.loc[(slice(None),slice(None),slice(None)),'dataset'] = 'daily mean'
daily_max_data.loc[(slice(None),slice(None),slice(None)),'dataset'] = 'daily max'
hourly_data.loc[(slice(None),slice(None),slice(None)),'dataset'] = 'hourly'

full_data = hourly_data.copy()
full_data = full_data.append(daily_mean_data)
full_data = full_data.append(daily_max_data)


## plotting data

full_data_no_emep = full_data.loc[(['random','startloss'],slice(None),slice(None)),]


sns.relplot(data=full_data_no_emep,x='slope',y='kendalltau_corr',row='dataset',col='chemical group',
    hue='chemical species',style='site_id',size='data loss scenario',kind='scatter')
#sns.relplot(data=full_data,x='slope',y='pearsonr_corr',row='stat',col='chemgroup',hue='spc',style='site_id',size='scenario',kind='scatter')

