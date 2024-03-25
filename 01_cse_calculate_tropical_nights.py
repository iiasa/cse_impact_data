# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 13:57:39 2022

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
# GCMs = ['MRI-ESM2-0']
# RCPs = ['ssp585']
# timestep = ''

# Set output directory
output_dir = ''

# Choose required tas and set variable
input_var = 'tasmin' # tas, tasmin, tasmax
output_var = 'tr20'
temp_threshold = 20

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_folder_with_repo\\required_files\\ISIMIP3b_GCM_GMT_1601_2100.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
year_range = 30
reference_period = [1974, 2004]


#%% Calculate tropical nights
#------------------------------------------------------------------------------

# Calculate min and max years for each GCM/RCP/threshold combination
years = cf.find_year_range(GCMs, RCPs, GMT_anomaly_file, GWLs, year_range)
    
for GCM, RCP in it.product(GCMs, RCPs):
    
    if len([val for key,val in years.items() if f'{GCM}_{RCP}' in key and len(val) != 0]) > 0:
        
        print(f'{GCM}_{RCP}')
        
        # Add one year at each end to allow calculation of tropical nights over a two year
        # period to avoid cutting off consecutive periods at the end of the year
        extended_years = [reference_period[0]-1, reference_period[1]+1]
        
        # Load data for reference period and convert to Celsius
        file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)
        tas_ref = cf.load_data(file_list, extended_years)
        tas_ref = cf.convert_to_Celsius(tas_ref)

        # Calculate indicator for all 31 years and average over all years
        tr20_ref = cfi.calculate_ctr20(tas_ref[input_var], temp_threshold, reference_period)
        tr20_ref_mean = tr20_ref.mean(dim='year')
        
        # Save annual data
        attributes = cfa.cust_attrs_tr20(GCM, RCP, 'historical', reference_period, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}_all')
        cf.write_output(tr20_ref, output_file, output_var, attributes)
       
        # Save averaged data
        attributes = cfa.cust_attrs_tr20(GCM, RCP, 'historical', reference_period, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}')
        cf.write_output(tr20_ref_mean, output_file, output_var, attributes)
     
    else:
        continue
    
    for gwl in GWLs:
        
        combination = f'{GCM}_{RCP}_{gwl}'
        
        # Check if the current threshold exists for the GCM/RCP combination
        if not years[combination]:
            continue
            
        else:
            
            # Add one year at each end to allow calculation of tropical nights over a two year
            # period to avoid cutting off consecutive periods at the end of the year
            extended_years = [years[combination][0]-1, years[combination][1]+1]
            
            # Load data
            file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)            
            tas = cf.load_data(file_list, extended_years)
            tas = cf.convert_to_Celsius(tas)
            
            # Calculate indicator for all 31 years and average over all years
            tr20 = cfi.calculate_ctr20(tas[input_var], temp_threshold, years[combination])
            tr20_mean = tr20.mean(dim='year')
            
            # Save annual data
            attributes = cfa.cust_attrs_tr20(GCM, RCP, gwl, years[combination], protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, 
                                                years[combination], f'{output_var}_all')
            cf.write_output(tr20, output_file, output_var, attributes)
            
            # Save averaged data
            attributes = cfa.cust_attrs_tr20(GCM, RCP, gwl, years[combination], protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, years[combination], f'{output_var}')
            cf.write_output(tr20_mean, output_file, output_var, attributes)

            