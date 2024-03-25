# -*- coding: utf-8 -*-
"""
Created on Fri May  6 09:02:58 2022

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
# GCMs = ['UKESM1-0-LL']
# RCPs = ['ssp126']
# timestep = ''

# Set output directory
output_dir = ''

# Choose required tas and set variable
input_var = 'pr' 
output_var = 'cdd'

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_folder_with_repo\\required_files\\ISIMIP3b_GCM_GMT_1601_2100.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
year_range = 30
reference_period = [1974, 2004]


#%% Calculate consecutive dry days 
#------------------------------------------------------------------------------

# Calculate min and max years for each GCM/RCP/threshold combination
years = cf.find_year_range(GCMs, RCPs, GMT_anomaly_file, GWLs, year_range)

for GCM, RCP in it.product(GCMs, RCPs):
       
    if len([val for key,val in years.items() if f'{GCM}_{RCP}' in key and len(val) != 0]) > 0:
    
        print(f'{GCM}_{RCP}')

        # Add one year at each end to allow calculation of cdd over a two year
        # period to avoid cutting off consecutive periods at the end of the year
        extended_years = [reference_period[0]-1, reference_period[1]+1]
        
        # Load data for reference period
        file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)
        pr_ref = cf.load_data(file_list, extended_years)
        
        # Calculate indicator for all 31 years and average over all years
        cdd_ref = cfi.calculate_cdd(pr_ref, reference_period)
        cdd_ref_mean = cdd_ref.mean(dim='year')
        
        # Save annual data
        attributes = cfa.cust_attrs_cdd(GCM, RCP, 'historical', reference_period, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}_all')
        cf.write_output(cdd_ref, output_file, output_var, attributes)
        
        # Save averaged data
        attributes = cfa.cust_attrs_cdd(GCM, RCP, 'historical', reference_period, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}')
        cf.write_output(cdd_ref_mean, output_file, output_var, attributes) 
               
    else:
        continue       
        
    for gwl in GWLs:
        
        combination = f'{GCM}_{RCP}_{gwl}'
        
        # Check if the current GWL exists for the GCM/RCP combination
        if not years[combination]:
            continue
            
        else:
            
            print(combination)
           
            # Add one year at each end to allow calculation of cdd over a two year
            # period to avoid cutting off consecutive periods at the end of the year
            extended_years = [years[combination][0]-1, years[combination][1]+1]

            # Load data
            file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)
            pr = cf.load_data(file_list, extended_years)
            
            # Calculate indicator for all 31 years and average over all years
            cdd = cfi.calculate_cdd(pr, years[combination])
            cdd_mean = cdd.mean(dim='year')
            
            # Save annual data
            attributes = cfa.cust_attrs_cdd(GCM, RCP, gwl, years[combination], protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, 
                                                years[combination], f'{output_var}_all')
            cf.write_output(cdd, output_file, output_var, attributes)
            
            # Save averaged data
            attributes = cfa.cust_attrs_cdd(GCM, RCP, gwl, years[combination], protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, 
                                                years[combination], f'{output_var}')
            cf.write_output(cdd_mean, output_file, output_var, attributes)     
        