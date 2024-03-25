# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 15:40:34 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import cse_functions as cf
import cse_functions_water as cfi
import cse_functions_attributes as cfa

#%% Settings
#------------------------------------------------------------------------------

# Set protocol
protocol = '3b' 
input_dir, GCMs, RCPs, timestep = cf.set_protocol(protocol)

# Manually overwrite GCMs/RCP/timestep here if required
# GCMs = ['UKESM1-0-LL']
# RCPs = ['ssp585']
# timestep = ''

# Set output directory
output_dir = ''

# Choose required tas and set variable
input_var = 'pr' 
output_var = 'precip'

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_folder_with_repo\\required_files\\ISIMIP3b_GCM_GMT_1601_2100.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
year_range = 30
reference_period = [1974, 2004]

#%% Calculate precipitation indicators
#------------------------------------------------------------------------------

# Calculate min and max years for each GCM/RCP/threshold combination
years = cf.find_year_range(GCMs, RCPs, GMT_anomaly_file, GWLs, year_range)
     
for GCM, RCP in it.product(GCMs, RCPs):
    
    if len([val for key,val in years.items() if f'{GCM}_{RCP}' in key and len(val) != 0]) > 0:
    
        print(f'{GCM}_{RCP}')
        
        # Load data for reference period
        file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)
        pr_ref = cf.load_data(file_list, reference_period)
        
        # Convert to mm and calculate quantiles
        pr_qnts_ref = cf.calculate_quantiles(pr_ref * 86400 , [0.95, 0.99], 50)
        
        # Calculate indicator for all 31 years and average over all years
        pr_ref = cfi.calculate_precipitation(pr_ref, pr_qnts_ref, reference_period)
        pr_ref_mean = pr_ref.mean(dim='year')
        
        # Save annual data
        attributes = cfa.cust_attrs_pr(GCM, RCP, 'historical', reference_period, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}_all')
        cf.write_output(pr_ref, output_file, output_var, attributes)
        
        # Save averaged data
        attributes = cfa.cust_attrs_pr(GCM, RCP, 'historical', reference_period, protocol)
        output_file = cf.create_output_file(output_dir, GCM, RCP, 'historical', 
                                            reference_period, f'{output_var}')
        cf.write_output(pr_ref_mean, output_file, output_var, attributes) 
        
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
            file_list = cf.create_file_list(input_dir, GCM, RCP, input_var, timestep)
            pr = cf.load_data(file_list, years[combination])

            # Calculate indicator for all 31 years and average over all years
            pr = cfi.calculate_precipitation(pr, pr_qnts_ref, years[combination])
            pr_mean = pr.mean(dim='year')
            
            # Save annual data
            attributes = cfa.cust_attrs_pr(GCM, RCP, gwl, years[combination], protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, 
                                                years[combination], f'{output_var}_all')
            cf.write_output(pr, output_file, output_var, attributes)
            
            # Save averaged data
            attributes = cfa.cust_attrs_pr(GCM, RCP, gwl, years[combination], protocol)
            output_file = cf.create_output_file(output_dir, GCM, RCP, gwl, 
                                                years[combination], f'{output_var}')
            cf.write_output(pr_mean, output_file, output_var, attributes)  