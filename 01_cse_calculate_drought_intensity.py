# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 14:29:56 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import cse_functions as cf
import cse_functions_water as cfi
import cse_functions_attributes as cfa
import xarray as xr

#%% Settings
#------------------------------------------------------------------------------

# Set protocol
protocol = '2b' 
input_dir, GCMs, RCPs, timestep = cf.set_protocol(protocol)
# Input dir needs to be set manually to directory with output data from GHMs
input_dir = 'path_to\\ISIMIP2b\\output'

# Manually overwrite GCMs/RCP/timestep here if required
# GCMs = ['MIROC5']
# RCPs = ['rcp60']
timestep = 'daily'

# Set output directory
output_dir = ''

# Specify GHMs and socs
GHMs = ['H08', 'LPJmL', 'MATSIRO']
soc = {'rcp26': ['rcp26soc'], 'rcp60': ['rcp60soc'], 'hist': ['histsoc']}

# Choose required tas and set variable
input_var = 'qtot'
output_var = 'dri_qtot'

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_folder_with_repo\\required_files\\ISIMIP2b_GCM_GMT_1661_2099.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
year_range = 30
reference_period = [1974, 2004]

# Load land area 
gridarea_path = 'path_to_folder_with_repo\\required_files\\gridarea05.nc'
land_area = xr.open_dataset(gridarea_path)

#%% Calculate drought intensity
#------------------------------------------------------------------------------

# Calculate min and max years for each GCM/RCP/threshold combination
years = cf.find_year_range(GCMs, RCPs, GMT_anomaly_file, GWLs, year_range)
     
for GHM, GCM, RCP in it.product(GHMs, GCMs, RCPs):
    
    if len([val for key,val in years.items() if f'{GCM}_{RCP}' in key and len(val) != 0]) > 0:
        
        print(f'{GHM}_{GCM}_{RCP}')
        
        # Load data for reference period
        file_list = cf.create_soc_file_list(input_dir, GCM, RCP, input_var, timestep, 'histsoc', 'histsoc', GHM)
        data_ref = cf.load_data(file_list, reference_period)        
        
        # If input variable is qtot, convert to correct unit
        if input_var == 'qtot':
            data_ref = data_ref * land_area['area'] * 1e-3
    
        # Calculate quantiles
        dri_q90_ref = cf.calculate_quantiles(data_ref, 0.9, 50)

        # Calculate indicator for all 31 years and average over all years
        dri = cfi.calculate_annual_drought_intensity(data_ref[input_var], 
                                                     dri_q90_ref, input_var)       
        dri_mean = dri.mean(dim='year')
        
        # Save annual data
        attributes = cfa.cust_attrs_di(GHM, GCM, RCP, 'histsoc', 'historical', reference_period, input_var, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'{output_var}_all', GHM)
        cf.write_output(dri, output_file, output_var, attributes)
         
        # Save averaged data
        attributes = cfa.cust_attrs_di(GHM, GCM, RCP, 'histsoc', 'historical', reference_period, input_var, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                output_var, GHM)
        cf.write_output(dri_mean, output_file, output_var, '')

        
    else:
        continue
    
    for gwl in GWLs:
    
        combination = f'{GCM}_{RCP}_{gwl}'
            
        # Check if the current threshold exists for the GCM/RCP combination
        if not years[combination]:
            continue
            
        else:
                          
            for SOC in soc[RCP]:
                
                print(f'{GHM}_{GCM}_{RCP}_{gwl}_{SOC}')
                
                # Load data
                file_list = cf.create_soc_file_list(input_dir, GCM, RCP, input_var, timestep, 'histsoc', SOC, GHM)
                data = cf.load_data(file_list, years[combination])
                
                # If input variable is qtot, convert to correct unit
                if input_var == 'qtot':
                    data = data * land_area['area'] * 1e-3
            
                # Calculate indicator for all 31 years and average over all years
                dri = cfi.calculate_annual_drought_intensity(data[input_var], dri_q90_ref, input_var)
                dri_mean = dri.mean(dim='year')

                # Save annual data                         
                attributes = cfa.cust_attrs_di(GHM, GCM, RCP, SOC, gwl, years[combination], input_var, protocol)
                output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 
                                                        SOC, gwl, years[combination],  
                                                        f'{output_var}_all', GHM)
                cf.write_output(dri, output_file, output_var, attributes)
                
                # Save averaged data
                attributes = cfa.cust_attrs_di(GHM, GCM, RCP, SOC, gwl, years[combination], input_var, protocol)
                output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 
                                                        SOC, gwl, years[combination], 
                                                        output_var, GHM)
                cf.write_output(dri_mean, output_file, output_var, attributes)