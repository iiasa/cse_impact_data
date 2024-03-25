# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:50:38 2022

@author: werning
"""

import sys
sys.path.append('path_to_folder_with_repo')
import itertools as it
import cse_functions as cf
import cse_functions_attributes as cfa
import xarray as xr
import os

#%% Settings
#------------------------------------------------------------------------------

# Set protocol
protocol = '2b' 
input_dir, GCMs, RCPs, timestep = cf.set_protocol(protocol)

# Input dir needs to be dir with raw wsi data
input_dir = ''

# Manually overwrite GCMs/RCP/timestep here if required
# GCMs = ['IPSL-CM5A-LR']
# RCPs = ['rcp26', 'rcp60']
timestep = 'daily'

# Set output directory
output_dir = ''

# Specify GHMs and socs
GHMs = ['H08', 'MATSIRO', 'LPJmL']
soc = {'rcp26': ['rcp26soc'], 'rcp60': ['rcp60soc'], 'hist': ['histsoc']}

# Set variable
input_var = 'wsi'

# Specify file with temperature anomalies, thresholds and year range
GMT_anomaly_file = 'path_to_folder_with_repo\\required_files\\ISIMIP2b_GCM_GMT_1661_2099.xlsx'
GWLs = [1.2, 1.5, 2.0, 2.5, 3.0, 3.5]
year_range = 30
reference_period = [1974, 2004]
supply_type = 'Inflow'  # this is default

# %% Calculate water stress index
# -----------------------------------------------------------------------------

# Calculate min and max years for each GCM/RCP/threshold combination
years = cf.find_year_range(GCMs, RCPs, GMT_anomaly_file, GWLs, year_range)

for GHM, GCM, RCP in it.product(GHMs, GCMs, RCPs):
    
    for SOC in soc[RCP]:
    
        print(f'{GHM}_{GCM}_{RCP}_{SOC}')
        
        # Open data period and reset dates
        data = xr.open_dataset(os.path.join(input_dir, f'{RCP}_{SOC}_Yearly', f'wsi_{GHM}_{GCM}_{RCP}_{SOC}_1971_2099.nc4'), decode_times=False)
        data['time'] =  [int(i)+1901 for i in data.time]
        data = data.rename({'time': 'year'})
        
        # Select required years and average over all years
        wsi_ref = data.sel(year=slice(reference_period[0], reference_period[1]))
        wsi_ref_mean = wsi_ref.mean(dim='year')
        
        # Save annual data
        attributes = cfa.cust_attrs_wsi(GHM, GCM, RCP, 'histsoc', input_var, 
                                        'historical', reference_period, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'{input_var}_all', GHM)
        cf.write_output(wsi_ref, output_file, input_var, attributes)
        
        # Save averaged data
        attributes = cfa.cust_attrs_wsi(GHM, GCM, RCP, 'histsoc', input_var, 
                                        'historical', reference_period, protocol)
        output_file = cf.create_soc_output_file(output_dir, GCM, RCP, 'histsoc', 
                                                'historical', reference_period, 
                                                f'{input_var}', GHM)
        cf.write_output(wsi_ref_mean, output_file, input_var, attributes)
        
        for gwl in GWLs:
            
            combination = f'{GCM}_{RCP}_{gwl}'
            
            if years[combination]:
                
                print(years[f'{GCM}_{RCP}_{gwl}'])
                
                # Select required years and average over all years
                wsi = data.sel(year=slice(years[combination][0], years[combination][1]))
                wsi_mean = wsi.mean(dim='year')
                
                # Save annual data 
                attributes = cfa.cust_attrs_wsi(GHM, GCM, RCP, SOC, input_var, 
                                                gwl, years[combination], protocol)
                output_file = cf.create_soc_output_file(output_dir, GCM, RCP, SOC, 
                                                        gwl, years[combination], 
                                                        f'{input_var}_all', GHM)
                cf.write_output(wsi, output_file, input_var, attributes)
                
                # Save averaged data
                attributes = cfa.cust_attrs_wsi(GHM, GCM, RCP, SOC, input_var, 
                                                gwl, years[combination], protocol)
                output_file = cf.create_soc_output_file(output_dir, GCM, RCP, SOC, 
                                                        gwl, years[combination], 
                                                        f'{input_var}', GHM)
                cf.write_output(wsi_mean, output_file, input_var, attributes)