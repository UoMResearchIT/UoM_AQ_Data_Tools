from selenium import webdriver
import time
import re


def process_pause(length):
	if(length=='short'):
		time.sleep(0.6)
	elif(length=='medium'):
		time.sleep(2.0)
	elif(length=='long'):
		time.sleep(4.0)
	else:
		time.sleep(5.0)


def search_and_select_radio_option(driver,want_data_type):
	radio_options = driver.find_elements_by_class_name('radio')
	while len(radio_options)>0:
		rad_opt = radio_options.pop()
		if(rad_opt.text==want_data_type):
			break
	rad_opt.find_element_by_class_name('checkround').click()


def find_number_of_results(driver):
	find_results = re.compile("(\d*) results found")
	count_str_list = find_results.findall(driver.page_source)
	return(int(count_str_list[0]))



def open_chrome_load_defra_website():
	driver = webdriver.Chrome()
	driver.get("https://uk-air.defra.gov.uk/data/data_selector")
	return(driver)



def select_hourly_networks(driver):
	driver.find_element_by_id('search-hourly-networks').click()
	driver.find_element_by_class_name('btn-start').click()
	process_pause('medium')  # wait for the page to load




def select_data_type(driver,want_data_type):

	driver.find_element_by_link_text("Select Data Type").click()
	process_pause('short')

	# get all the options in the list of radio buttons
	search_and_select_radio_option(driver,want_data_type)
	process_pause('short')

	# search for the data type we are interested in, and select it
	search_and_select_radio_option(driver,want_data_type)
	process_pause('short')

	driver.find_element_by_id('submit').click()
	process_pause('medium')



def select_date_range(driver,start_date,end_date):
	want_data_type = 'Custom Date (below)'
	driver.find_element_by_link_text('Select Date Range').click()
	process_pause('short')

	# search for the custom date radio button, and click it
	search_and_select_radio_option(driver,want_data_type)
	process_pause('short')

	date_from = driver.find_element_by_xpath("//*[@id='f_date_started']")
	date_from.clear()
	date_from.send_keys(start_date)
	process_pause('short')

	date_to = driver.find_element_by_xpath("//*[@id='f_date_ended']")
	date_to.clear()
	date_to.send_keys(end_date)
	process_pause('short')

	driver.find_element_by_id('submit').click()
	process_pause('medium')


def select_monitoring_site(driver,work_startend):
	driver.find_element_by_link_text('Select Monitoring Sites').click()
	process_pause('short')
	driver.find_element_by_link_text('Monitoring Network').click()
	process_pause('short')

	site_menu = driver.find_element_by_name('f_group_id')
	if(site_menu.get_attribute('value') != '4'):
		print('wrong network!')
	else:
		print('AURN selected')
	process_pause('short')

	driver.find_element_by_id('submit').click()
	process_pause('short')

	# pull out information on number of sites, and set range to span
	nsites = find_number_of_results(driver)
	full_range = list(range(1,nsites))

	# determine working range, limiting it to less than full range
	work_range = list(range(min(work_startend[0],nsites),min(work_startend[1],nsites)))

	# make sure all options are deselected
	for val in full_range:
		test_opt = driver.find_element_by_xpath("//*[@id='f_site_id']/option["+str(val)+"]")
		if(test_opt.is_selected()):
			test_opt.click()
	

	# selection the options we want
	for val in work_range:
		test_opt = driver.find_element_by_xpath("//*[@id='f_site_id']/option["+str(val)+"]")
		if(not test_opt.is_selected()):
			test_opt.click()
	process_pause('medium')

	driver.find_element_by_id('submit').click()
	process_pause('medium')



def select_pollutants(driver,target_polls):
	driver.find_element_by_link_text('Select Pollutants').click()
	process_pause('short')
	driver.find_element_by_link_text('Monitoring Network').click()
	process_pause('short')

	site_menu = driver.find_element_by_name('f_network_id[]')
	if(site_menu.get_attribute('value') != '4'):
		print('wrong network!')
	else:
		print('AURN selected')
	process_pause('short')

	driver.find_element_by_id('submit').click()
	process_pause('short')

	# make sure all options are deselected, then reselect pollutants of interest
	val = 1
	while val > 0:
		try:
			test_opt = driver.find_element_by_xpath("//*[@id='f_parameter_id']/option["+str(val)+"]")
			if(test_opt.is_selected()):
				test_opt.click()
			if(test_opt.get_attribute('value') in target_polls):
				test_opt.click()
			val += 1
		except:
			val = -1
	process_pause('short')

	driver.find_element_by_id('submit').click()
	process_pause('medium')


def select_output_type(driver,email_string):

	driver.find_element_by_link_text('Selected Output Type').click()
	process_pause('short')

	output_select = driver.find_element_by_name('output_type_selector')
	opt_list = output_select.find_element_by_class_name('radio').find_elements_by_class_name('radio')
	while len(opt_list)>0:
		output_opt = opt_list.pop()
		if(output_opt.text=='Data to Email Address (CSV)'):
			break

	output_opt.find_element_by_class_name('checkround').click()
	process_pause('short')

	driver.find_element_by_id('f_email').send_keys(email_string)
	process_pause('short')

	driver.find_element_by_xpath('//*[@id="gdpr_tc"]/label/span').click()
	process_pause('short')

	driver.find_element_by_id('submit').click()
	process_pause('medium')


def submit_final_form(driver):
	process_pause('long')
	driver.find_element_by_id('getDataButton').click()
	process_pause('medium')

def reset_search_form(driver):
	driver.find_element_by_class_name('btn-sm').click()
	process_pause('medium')

def run_workflow(driver,process_list,process_start_index,process_end_index,start_date,end_date,email_string):
	
	process_range = list(range(process_start_index,min(process_end_index,len(process_list))))
	
	for val in process_range:
		data_type = process_list[val]['data_type']
		sites_range = process_list[val]['sites_range']
		pollutants = process_list[val]['pollutants']
		print(f'Extraction position {val}: {data_type} of {pollutants} for sites {sites_range}')
		select_hourly_networks(driver)
		select_data_type(driver,data_type)
		select_date_range(driver,start_date,end_date)
		select_monitoring_site(driver,sites_range)
		select_pollutants(driver,pollutants)
		select_output_type(driver,email_string)
		
		submit_final_form(driver)
		reset_search_form(driver)


def build_workflow_process_list(data_types,sites_list,pollutants_list):

	process_list = []
	for data_type in data_types:
		for sites_range in sites_list:
			for pollutants in pollutants_list:
				process_list.append({"data_type":data_type, "sites_range":sites_range, "pollutants":pollutants})

	return(process_list)


if __name__ == '__main__':

	##### start of settings section

	email_string='ENTER YOUR EMAIL ADDRESS HERE'
	start_date = '01/01/2019'
	end_date   = '31/12/2019'
	
	
	# data types to extract - select daily data types only
	data_types = ['Daily Max','Daily Mean']

	# Numerical site numbers - each pair of numbers should be start & end points (for 
	# breaking the list into chunks).
	#
	# The DEFRA website limits the size of data arrays which can be downloaded at one time.
	# Trail and error has shown 90 sites to be a reasonable number (when a single pollutant
	# is selected).
	#
	# The starting index is 2 (1st index is "all", do *not* use this!)
	# The final end point should be large enough to cover all possible sites (the 
	# minimum of this and the actual list size will be used by the script)
	sites_list = []
	sites_list.append([2,90])
	sites_list.append([91,200])

	# List the (grouped) pollutants required. Each group can contain more than one 
	# pollutant if desired, but if doing so then be careful the size of the site list, to
	# avoid requesting too large a data array). 
	pollutants_list = []
	pollutants_list.append(['NO2'])       # NO2
	pollutants_list.append(['NOXasNO2'])  # NOx, reported as mass of NO2
	pollutants_list.append(['O3'])        # O3
	pollutants_list.append(['SO2'])       # SO2
	pollutants_list.append(['GE10'])      # PM10
	pollutants_list.append(['PM25'])      # PM2.5
	
	# Use a start index > 0 for restarting a failed process from the first failed process
	#
	# The end index should be greater than [# data_types] x [# sites_list] x [# pollutants list]
	# unless you need to select a subset of the list.
	process_start_index = 0
	process_end_index = 200

	##### end of settings section

	# build processing list
	process_list = build_workflow_process_list(data_types,sites_list,pollutants_list)
	
	# echo the settings within the process list to be used
	process_range = list(range(process_start_index,min(process_end_index,len(process_list))))
	for val in process_range:
		print(val, process_list[val])

	# load the DEFRA website
	driver = open_chrome_load_defra_website()
	process_pause('long')

	# Run the workflow
	run_workflow(driver,process_list,process_start_index,process_end_index,start_date,end_date,email_string)
