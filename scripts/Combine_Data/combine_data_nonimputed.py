import pandas as pd


aurn_file='../AURN_Data_Download/AURN_data/aurn_processed_daily_2016-2019.csv'
emep_file='../EMEP_Data_Extraction/EMEP_data/emep_daily_data_2016-2019.csv'
poll_file='../MEDMI_Data_Download/full_data/pollen_2016-2019.csv'
met_file='../Data_Processing/MEDMI_Met_data/Met_ppd_daily_mean_max_temp_RH_pres_2016-2019_no_imputation.csv'

outfile='Combined_dataset/turing_aq_daily_met_pollen_pollution_data.csv'


aurn_data = pd.read_csv(aurn_file,index_col=['timestamp','site_id'])
emep_data = pd.read_csv(emep_file,index_col=['timestamp','site_id'])
poll_data = pd.read_csv(poll_file,index_col=['timestamp','site_id'])
met_data = pd.read_csv(met_file,index_col=['timestamp','site_id'])


combined_data = aurn_data.copy()
combined_data = combined_data.merge(emep_data, how='outer', left_index=True, right_index=True)
combined_data = combined_data.merge(poll_data, how='outer', left_index=True, right_index=True)
combined_data = combined_data.merge(met_data, how='outer', left_index=True, right_index=True)

combined_data.to_csv(outfile,index=True,header=True,float_format='%.2f')
