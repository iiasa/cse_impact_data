# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 10:23:08 2021

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import cse_functions as cf
import cse_functions_indicators as cfi
import cse_functions_attributes as cfa

#%% Settings
#------------------------------------------------------------------------------

# Set protocol
protocol = '3b' 
input_dir, GCMs, RCPs, timestep = cf.set_protocol(protocol)

# Manually overwrite GCMs/RCP/timestep here if required
# GCMs = ['MIROC5']
# RCPs = ['rcp26', 'rcp60', 'rcp85']
# timestep = ''

# Set output directory
output_dir = ''
  
# Choose required tas and set variable
input_var = 'tas' # tas, tasmin, tasmax
output_var = 'sdd'

# Set balance temperatures
balance_temperature_heating_low = 18.0
balance_temperature_cooling_low = 24.0

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_folder_with_repo\\required_files\\ISIMIP3b_GCM_GMT_1601_2100.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
year_range = 30
reference_period = [1974, 2004]

#%% Calculate simple degree days temperature
#------------------------------------------------------------------------------

# Calculate min and max years for each GCM/RCP/threshold combination
years = cf.find_year_range(GCMs, RCPs, GMT_anomaly_file, GWLs, year_range)
     
for GCM, RCP in it.product(GCMs, RCPs):
    
    if len([val for key,val in years.items() if f'{GCM}_{RCP}' in key and len(val) != 0]) > 0:
    
        print(f'{GCM}_{RCP}')
        
        # Load data for reference period and convert to Celsius
        file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)
        tas_ref = cf.load_data(file_list, reference_period)
        tas_ref = cf.convert_to_Celsius(tas_ref)      
        
        # Calculate indicator for all 31 years and average over all years
        simple_degree_days_ref = cfi.calculate_sdd(tas_ref, balance_temperature_cooling_low, 
                                               balance_temperature_heating_low)
        simple_degree_days_ref_mean = simple_degree_days_ref.mean(dim='year') 
        
        # Save annual data
        attributes = cfa.cust_attrs_sdd(GCM, RCP, 'historical', reference_period, 
                                        balance_temperature_cooling_low, 
                                        balance_temperature_heating_low, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, 
                                            f'{output_var}{input_var.replace("tas", "")}_{str(balance_temperature_cooling_low).replace(".", "p")}_all')
        cf.write_output(simple_degree_days_ref, output_file, output_var, attributes)
        
        # Save averaged data
        attributes = cfa.cust_attrs_sdd(GCM, RCP, 'historical', reference_period, 
                                        balance_temperature_cooling_low, 
                                        balance_temperature_heating_low, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, 
                                            f'{output_var}{input_var.replace("tas", "")}_{str(balance_temperature_cooling_low).replace(".", "p")}')
        cf.write_output(simple_degree_days_ref_mean, output_file, output_var, attributes)
        
    else:
        continue
    
    for gwl in GWLs:
        
        combination = f'{GCM}_{RCP}_{gwl}'
        
        # Check if the current threshold exists for the GCM/RCP combination
        if not years[combination]:
            continue
            
        else:
            
            print(combination)
           
            # Load data
            tas = cf.load_data(file_list, years[combination])
            tas = cf.convert_to_Celsius(tas)      
            
            # Calculate indicator for all 31 years and average over all years
            simple_degree_days = cfi.calculate_sdd(tas, balance_temperature_cooling_low, 
                                                   balance_temperature_heating_low)
            simple_degree_days_mean = simple_degree_days.mean(dim='year') 
                        
            # Save annual data
            attributes = cfa.cust_attrs_sdd(GCM, RCP, gwl, years[combination], 
                                            balance_temperature_cooling_low, 
                                            balance_temperature_heating_low, protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl,
                                                years[combination], 
                                                f'{output_var}{input_var.replace("tas", "")}_{str(balance_temperature_cooling_low).replace(".", "p")}_all')
            cf.write_output(simple_degree_days, output_file, output_var, attributes)
            
            # Save averaged data
            attributes = cfa.cust_attrs_sdd(GCM, RCP, gwl, years[combination], 
                                            balance_temperature_cooling_low, 
                                            balance_temperature_heating_low, protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, 
                                                years[combination], 
                                                f'{output_var}{input_var.replace("tas", "")}_{str(balance_temperature_cooling_low).replace(".", "p")}')
            cf.write_output(simple_degree_days_mean, output_file, output_var, attributes)  