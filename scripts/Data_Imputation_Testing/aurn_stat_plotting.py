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

#chemical_scenarios = ['NO2','NOXasNO2','O3','PM10','PM2.5','SO2','combined','combined_and_so2']
chemical_scenarios = ['NO2','NOXasNO2','O3','PM10','PM2.5','SO2','combined_and_so2']
chemical_dictionary = {'NO2':'individual','NOXasNO2':'individual','O3':'individual','PM10':'individual',
                       'PM2.5':'individual','SO2':'individual','combined':'combined','combined_and_so2':'combined_and_so2'}
chemical_mixes = ['individual','combined','combined_and_so2']

chemical_dictionary = {'NO2':'individual','NOXasNO2':'individual','O3':'individual','PM10':'individual',
                       'PM2.5':'individual','SO2':'individual','combined_and_so2':'combined'}
chemical_mixes = ['individual','combined']

sample_scenarios = ['random','random_with_EMEP','startloss','startloss_with_EMEP','startloss0.25','startloss0.75',
                    'startloss0.25_with_EMEP','startloss0.75_with_EMEP']

sample_scenarios = ['random','random_with_EMEP','startloss','startloss_with_EMEP']


stats_list = ['kendalltau_corr','pearsonr_corr','spearmanr_corr','slope']
index_list = ['site_id','spc']
col_hourly_list = stats_list + index_list
col_daily_list = col_hourly_list + ['stat']


chemical_pairs = {'NO2':'NOx','NOXasNO2':'NOx','O3':'SO2 O3','PM10':'PM','PM2.5':'PM','SO2':'SO2 O3'}

column_substitutes = {'spc':'chemical species','spearmanr_corr':"Spearman's rank correlation",'slope':'slope of fit'}




def read_datafile(full_data,filepath,index_list_internal,col_list_internal,scenario_string,chemical_string,chem_pairs_internal):
    try:
        temp_data = pd.read_csv(filepath,usecols=col_list_internal)
        temp_data['data loss scenario'] = scenario_string
        if scenario_string.split('_')[-1] == 'EMEP':
            chem_string_wrap = '{}_EMEP'
        else:
            chem_string_wrap = '{}'
        temp_data['chemical mix'] = chem_string_wrap.format(chemical_string)
        temp_data['chemical group'] = [chem_pairs_internal[chemkey] for chemkey in temp_data['spc']]
        temp_data = temp_data.rename(columns=column_substitutes)
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


### isolate particular scenarios

full_data_random = full_data.loc[(['random','random_with_EMEP'],slice(None),slice(None)),]
full_data_startloss = full_data.loc[(['startloss','startloss_with_EMEP'],slice(None),slice(None)),]

full_data_startloss_only = full_data.loc[(['startloss_with_EMEP','startloss0.25_with_EMEP','startloss0.75_with_EMEP'],'combined_and_so2',slice(None)),]
full_data_random_only = full_data.loc[(['random'],slice(None),slice(None)),]


paper_data_random = full_data.loc[(['random','random_with_EMEP'],slice(None),['NO2','PM10']),]

paper_data_random = full_data.loc[(['random','random_with_EMEP'],['individual','combined','combined_EMEP'],['NO2','PM10']),]

paper_data_startloss = full_data.loc[(['startloss','startloss_with_EMEP'],['individual','combined','combined_EMEP'],['NO2','PM10']),]

## plotting data

sns.set_style('ticks')
sns.set_context('talk',font_scale=0.9)
sns.despine()

col_order = ['NOx','SO2 O3','PM']
row_order = ['hourly','daily mean','daily max']
style_order = ['MY1','PT4','LEED','MID']

gplot = sns.relplot(data=full_data_random,x='slope of fit',y="Spearman's rank correlation",row='dataset',col='chemical group',
    row_order = row_order, col_order = col_order, style_order = style_order, sizes=(40,100),
    hue='chemical species',style='site_id',size='data loss scenario',kind='scatter',legend='full',height=3,aspect=1.33)
gplot.set(xlim=(-0.01,1.2),ylim=(-0.01,1.01))
#sns.relplot(data=full_data,x='slope',y='pearsonr_corr',row='stat',col='chemgroup',hue='spc',style='site_id',size='scenario',kind='scatter')

gplot = sns.relplot(data=full_data_random_only,x='slope of fit',y="Spearman's rank correlation",row='dataset',col='chemical group',
    row_order = row_order, col_order = col_order, style_order = style_order, sizes=(40,100),
    hue='chemical species',style='site_id',size='chemical mix',kind='scatter',legend='full',height=3,aspect=1.33)
gplot.set(xlim=(-0.01,1.2),ylim=(-0.01,1.01))


## 
size_order = ['startloss0.25_with_EMEP','startloss_with_EMEP','startloss0.75_with_EMEP']
gplot = sns.relplot(data=full_data_startloss_only,x='slope of fit',y="Spearman's rank correlation",row='dataset',col='chemical group',
    row_order = row_order, col_order = col_order, style_order = style_order, size_order = size_order, sizes=(40,100),
    hue='chemical species',style='site_id',size='data loss scenario',kind='scatter',legend='full',height=3,aspect=1.33)
gplot.set(xlim=(-0.01,1.1),ylim=(-0.01,1.01))



### paper
gplot = sns.relplot(data=paper_data_startloss[paper_data_startloss['dataset']=='hourly'],x='slope of fit',y="Spearman's rank correlation",col='chemical species',
    hue='chemical mix',style='site_id',kind='scatter',legend='full',height=3,aspect=1.33)
gplot.set(xlim=(-0.01,1.1),ylim=(-0.01,1.01))

