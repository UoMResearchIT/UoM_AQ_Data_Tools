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

daily_file = 'met_daily_correlation_stats.csv'
hourly_file = 'met_hourly_correlation_stats.csv'

path_head = 'met_stats/test_stats_{}_{}'

met_scenarios = ['sitesA','sitesB','sitesC']
chemical_dictionary = {'NO2':'individual','NOXasNO2':'individual','O3':'individual','PM10':'individual',
                       'PM2.5':'individual','SO2':'individual','combined':'combined','combined_and_so2':'combined_and_so2'}
chemical_mixes = ['individual','combined','combined_and_so2']

sample_scenarios = ['random_50','random_25','random_75','startloss_25','startloss_50','startloss_75','endloss_25','endloss_50','endloss_75']
sample_scenarios = ['random_50','random_25','startloss_25','startloss_50','endloss_25','endloss_50']

stats_list = ['kendalltau_corr','pearsonr_corr','spearmanr_corr','slope']
index_list = ['site_id','spc']
col_hourly_list = stats_list + index_list
col_daily_list = col_hourly_list + ['stat']

column_substitutes = {'spc':'met variable','relativehumidity':'relative humidity',
                      'spearmanr_corr':"Spearman's rank correlation",'slope':'slope of fit'}




def read_datafile(full_data,filepath,index_list_internal,col_list_internal,scenario_string):
    try:
        temp_data = pd.read_csv(filepath,usecols=col_list_internal)
        temp_data['data loss scenario'] = scenario_string
        temp_data = temp_data.rename(columns=column_substitutes)
        temp_data['site_id'] = [str(siteid) for siteid in temp_data['site_id']]
        return(full_data.append(temp_data))
    except:
        print('could not process file {}'.format(filepath))
        return(full_data)


daily_data  = pd.DataFrame()
hourly_data = pd.DataFrame()


for scenario in sample_scenarios:
    for metscen in met_scenarios:
        base_dir = Path(path_head.format(metscen,scenario))

        if base_dir.is_dir():
            hourly_path = base_dir.joinpath(hourly_file) 
            daily_path  = base_dir.joinpath(daily_file)
            daily_data = read_datafile(daily_data,daily_path,index_list,col_daily_list,scenario)
            hourly_data = read_datafile(hourly_data,hourly_path,index_list,col_hourly_list,scenario)
        else:
            print('{} directory does not exist'.format(base_dir))

daily_mean_data = daily_data[daily_data['stat']=='mean']
daily_max_data = daily_data[daily_data['stat']=='max']
daily_mean_data.set_index(['data loss scenario','met variable'],inplace=True)
daily_max_data.set_index(['data loss scenario','met variable'],inplace=True)
hourly_data.set_index(['data loss scenario','met variable'],inplace=True)

daily_mean_data.loc[(slice(None),slice(None)),'dataset'] = 'daily mean'
daily_max_data.loc[(slice(None),slice(None)),'dataset'] = 'daily max'
hourly_data.loc[(slice(None),slice(None)),'dataset'] = 'hourly'

full_data = hourly_data.copy()
full_data = full_data.append(daily_mean_data)
full_data = full_data.append(daily_max_data)


## plotting data

sns.set_style('ticks')
sns.set_context('talk',font_scale=0.9)
sns.despine()


gplot = sns.relplot(data=full_data,x='slope of fit',y="Spearman's rank correlation",row='dataset',col='met variable',
    hue='site_id',style='data loss scenario',kind='scatter',height=3,aspect=1.33)
gplot.set(xlim=(-0.01,1.2),ylim=(-0.01,1.01))


gplot = sns.relplot(data=full_data.loc[(slice(None),'relativehumidity'),],x='slope of fit',y="Spearman's rank correlation",col='dataset',
    hue='site_id',style='data loss scenario',kind='scatter',height=3,aspect=1.33)
gplot.set(xlim=(0.39,1.01),ylim=(0.39,1.01))
