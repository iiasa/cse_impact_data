# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 09:24:53 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import cse_functions as cf
import cse_functions_indicators as cfi
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
# GCMs = ['GFDL-ESM2M']
# RCPs = ['rcp26', 'rcp60']
# RCPs = ['rcp26']
timestep = 'daily'

# Set output directory
output_dir = ''

# Specify GHMs and socs
GHMs = ['H08', 'LPJmL', 'MATSIRO']
soc = {'rcp26': ['rcp26soc'], 'rcp60': ['rcp60soc'], 'hist': ['histsoc']}

# Set variable
input_var = 'dis'

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_folder_with_repo\\required_files\\ISIMIP2b_GCM_GMT_1661_2099.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
year_range = 30
reference_period = [1974, 2004]

# Load land area
gridarea_path = 'path_to_folder_with_repo\\required_files\\gridarea05.nc'
land_area = xr.open_dataset(gridarea_path)

# %% Calculate seasonality and inter-annual variability
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
        
        # Calculate indicator for all 31 years and average over all years
        seas_ref, iavar_ref, monthly_mean_ref = cfi.calculate_seasonality(data_ref[input_var])
        seas_annual_ref, iavar_annual_ref = cfi.calculate_annual_seasonality(data_ref[input_var])

        # Save annual data
        attributes = cfa.cust_attrs_seasonality(GHM, GCM, RCP, 'histsoc', input_var, 
                                                'historical', reference_period, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'seasonality_{input_var}_all', GHM)
        cf.write_output(seas_annual_ref, output_file, 'seasonality', '')
    
        attributes = cfa.cust_attrs_iavar(GHM, GCM, RCP, 'histsoc', input_var, 
                                          'historical', reference_period, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'ia_var_{input_var}_all', GHM)
        cf.write_output(iavar_annual_ref, output_file, 'ia_var', '')
         
        # Save averaged data
        attributes = cfa.cust_attrs_seasonality(GHM, GCM, RCP, 'histsoc', input_var, 
                                                'historical', reference_period, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'seasonality_{input_var}', GHM)
        cf.write_output(seas_ref, output_file, 'seasonality', attributes)
    
        attributes = cfa.cust_attrs_iavar(GHM, GCM, RCP, 'histsoc', input_var, 
                                          'historical', reference_period, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'ia_var_{input_var}', GHM)
        cf.write_output(iavar_ref, output_file, 'ia_var', attributes)
    
        attributes = cfa.cust_attrs_monthly_mean(GHM, GCM, RCP, 'histsoc', input_var, 
                                                 'historical', reference_period, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'monthly_mean_{input_var}', GHM)
        cf.write_output(monthly_mean_ref, output_file, 'monthly_mean', attributes)      
        
        
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
                file_list = cf.create_soc_file_list(input_dir, GCM, RCP, input_var, 
                                                    timestep, 'histsoc', SOC, GHM)
                data = cf.load_data(file_list, years[combination])
                
                # If input variable is qtot, convert to correct unit
                if input_var == 'qtot':
                    data = data * land_area['area'] * 1e-3
    
                # Calculate indicator 
                seas, iavar, monthly_mean = cfi.calculate_seasonality(data[input_var])
    
                # Save data
                attributes = cfa.cust_attrs_seasonality(GHM, GCM, RCP, SOC, 
                                                        input_var, gwl, years[combination], 
                                                        protocol)
                output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 
                                                        SOC, gwl, years[combination], 
                                                        f'seasonality_{input_var}', GHM)
                cf.write_output(seas, output_file, 'seasonality', attributes)
            
                attributes = cfa.cust_attrs_iavar(GHM, GCM, RCP, SOC, input_var, 
                                                  gwl, years[combination], protocol)
                output_file = cf.create_soc_output_file(output_dir, GCM, RCP, SOC, 
                                                        gwl, years[combination], 
                                                        f'ia_var_{input_var}', GHM)
                cf.write_output(iavar, output_file, 'ia_var', attributes)
            
                attributes = cfa.cust_attrs_monthly_mean(GHM, GCM, RCP, SOC, input_var, 
                                                         gwl, years[combination], protocol)
                output_file = cf.create_soc_output_file(output_dir, GCM, RCP, SOC, 
                                                        gwl, years[combination], 
                                                        f'monthly_mean_{input_var}', GHM)
                cf.write_output(monthly_mean, output_file, 'monthly_mean', attributes)
